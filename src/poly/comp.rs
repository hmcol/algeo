use std::marker::PhantomData;

use itertools::{EitherOrBoth, Itertools};

use crate::core::num::Field;

use super::mdeg::MultiDegree;
use super::ord::MonomialOrder;
use super::elts::*;

/// Handler for computations related to finding Gröbner bases.
///
/// The field `F` is the field of the polynomials it works over and the monomial order `O` is what it uses to determine leading terms
pub struct Computer<F: Field, O: MonomialOrder> {
    _marker: PhantomData<(F, O)>,
}

impl<F: Field, O: MonomialOrder> Computer<F, O> {
    /// Returns the polynomial with the same terms as `f`, sorted by the monomial order `Self::O`.
    pub fn sort_terms(f: &Polynomial<F>) -> Polynomial<F> {
        // clean operation:
        // trust that `f` has clean terms, so reordering will not spoil
        Polynomial::new_unchecked(
            f.terms()
                .sorted_by(|s, t| O::cmp(&s.mdeg, &t.mdeg))
                .cloned()
                .collect(),
        )
    }

    /// Returns the leading term of `f`, if `f` is nonzero
    ///
    /// Returns `None` if `f` is zero (has no terms), corresponding to the convention $LT(0) = 0$.
    ///
    /// Ideally, should also return none if it only has zero terms; but also we would ideally never have zero terms, so go figure.
    pub fn leading_term(f: &Polynomial<F>) -> Option<&Term<F>> {
        f.terms().max_by(|s, t| O::cmp(&s.mdeg, &t.mdeg))
    }

    /// Returns the coefficient of the leading term of `f`, if `f` is nonzero
    ///
    /// Returns `None` if `f` is zero (has no terms), corresponding to the convention $LT(0) = 0$.
    pub fn leading_coef(f: &Polynomial<F>) -> Option<F> {
        Some(Self::leading_term(f)?.coef)
    }

    /// Performs general polynomial division on `f` with the divisors `divs`; returns the pair `(remainder, quotients)`.
    ///
    /// Let `f` = $f$ and `divs` = $\\{g_1, \dots, g_m\\}$.
    ///
    /// Returns $(r, \\{q_1, \dots, q_m\\})$, where
    /// $$ f = q_1g_1 + \cdots + q_mg_m + r, $$
    /// and the multidegree of each $q_ig_i$ is no more than te multidegree of $f$.
    ///
    /// If it happens that $r = 0$, then we can conclude $f \in I = \langle g_1, \dots, g_m \rangle$.
    ///
    /// If the divisors form a Gröbner basis for $I$, then the order of the $g_i$'s does not affect the result. In particular $r = 0$ if and only if $f \in I$.
    ///
    /// If the divisors do no form a Gröbner basis, this is not necessarily the case.
    pub fn divide(
        f: &Polynomial<F>,
        divs: &[Polynomial<F>],
    ) -> (Polynomial<F>, Vec<Polynomial<F>>) {
        let m = divs.len();
        let mut quotients = vec![Polynomial::<F>::zero(); m];
        let mut remainder = Polynomial::<F>::zero();

        let mut f = f.clone();
        // let divs: Vec<Polynomial<F>> = divs.iter().filter(|g| !g.is_zero()).cloned().collect();

        'outer: while let Some(lt_f) = Self::leading_term(&f).cloned() {
            // f still has (nonzero) terms

            for (g, q) in divs.iter().zip(quotients.iter_mut()) {
                if let Some(lt_g) = Self::leading_term(&g) {
                    if let Some(a) = lt_f.try_div(lt_g) {
                        // case 1: `lt_f` is divisible by `lt_g`
                        // - `a` is the term such that `lt_f == a * lt_g`

                        // add `a` to `g`'s quotient `q`
                        *q = &*q + &a;

                        // eliminate `lt_f == a * lt_g` from `f`
                        f = f - a * g;

                        continue 'outer;
                    }
                }
            }

            // case 2: LT(f) is not divisible by any LT(g_1), ..., LT(g_m)
            // - subtract the leading term of f and add it to the remainder

            // sort the terms of f by the current order
            f.terms.sort_by(|s, t| O::cmp(&s.mdeg, &t.mdeg));

            // ofter sorting, last term of f is the leading term, LT(f)
            if let Some(lt_f) = f.terms.pop() {
                // first

                remainder = remainder + lt_f;
            }
        }

        // return
        (remainder, quotients)
    }

    /// Returns the monic least common multiple of the given terms
    ///
    /// Given $s = c x_1^{a_1} \dots x_n^{a_n}$ and $t = d x_1^{b_1} \dots x_n^{b_n}$, then `monic_lcm(s, t)` returns
    ///
    /// If $\partial s = (a_1, \dots, a_n)$ and $\partial t = (b_1, \dots, b_n)$, then
    ///
    /// $$ \text{monic-lcm}(s, t) = x_1^{\max(a_1, b_1)} \cdots x_n^{\max(a_n, b_n)}. $$
    ///
    /// won't work/doesn't make sense for for negative degrees
    pub fn monic_lcm(s: &Term<F>, t: &Term<F>) -> Term<F> {
        if s.is_zero() || t.is_zero() {
            // Not 100% confident this is the correct behavior
            // may want to return opposite value when one is zero
            return Term::zero();
        }
        // `s` and `t` guaranteed to be nonzero

        // clean operation:
        // - trust that `s.mdeg` and `t.mdeg` are clean
        // - taking either or max will not spoil
        Term::monic(MultiDegree(
            s.mdeg
                .degs()
                .zip_longest(t.mdeg.degs())
                .map(|pair| match pair {
                    EitherOrBoth::Both(a, b) => *a.max(b),
                    EitherOrBoth::Left(a) => *a,
                    EitherOrBoth::Right(b) => *b,
                })
                .collect(),
        ))
    }

    /// Returns a combination of `f` and `g` with their leading terms cancelled.
    ///
    /// In particular, if $M$ is the monic lcm of $\operatorname{LT}(f)$ and $\operatorname{LT}(g)$, `reduce(f, g)` returns
    ///
    /// $$ S(f, g) = \frac{M}{\operatorname{LT}(f)} f - \frac{M}{\operatorname{LT}(g)} g. $$
    pub fn try_reduce(f: &Polynomial<F>, g: &Polynomial<F>) -> Option<Polynomial<F>> {
        let lt_f = Self::leading_term(f)?;
        let lt_g = Self::leading_term(g)?;
        let lcm = Self::monic_lcm(&lt_f, &lt_g);

        Some(lcm.try_div(&lt_f)? * f - lcm.try_div(&lt_g)? * g)
    }

    /// Test whether `generators` is a Gröbner basis using Buchberger's criterion
    ///
    /// `generators` = $G = \{g_1, \dots, g_m\}$
    ///
    /// Buchberger's criterion is as follows: If $I = \langle G \rangle$ is a nonzero ideal, then $G$ is a Gröbner basis for $I$ if and only if $S(G_i, g_j) \equiv 0 \mod G$ for all $1 \leq i, j \leq m$
    pub fn buchberger_criterion(generators: &[Polynomial<F>]) -> bool {
        let m = generators.len();

        for (i, j) in (0..m).cartesian_product(0..m).filter(|&(f, g)| f != g) {
            if let Some(h) = Self::try_reduce(&generators[i], &generators[j]) {
                // h = S(g_i, g_j)

                if !Self::divide(&h, generators).0.is_zero() {
                    // h != 0 (mod G)

                    return false;
                }
            } else {
                panic!("should not happen. yeah i know, joseph, 'parse don't validate'");
            }
        }

        true
    }

    #[allow(clippy::mut_range_bound)]
    pub fn buchberger_algorithm(generators: &mut Vec<Polynomial<F>>) {
        let mut m = generators.len();
        let mut i = 0;

        'extend: while i < m {
            for j in 0..i {
                // remainder of `S(g_i, g_j)` after division by `generator`
                let (r, _) = Self::divide(
                    &Self::try_reduce(&generators[i], &generators[j])
                        .expect("should not happen. yeah i know, joseph, 'parse don't validate'"),
                    generators,
                );

                if let Some(lc_r) = Self::leading_coef(&r) {
                    // `r` is nonzero with leading coef `lc_r`
                    // in other words, S(g_i, g_j) = 0 mod G

                    // add remainder to the generators
                    generators.push(&r / lc_r);
                    m += 1;

                    // start the process over
                    
                    i = 0;
                    continue 'extend;
                }
            }
            i += 1;
        }

        #[inline]
        fn pairs_iter(m: usize) -> impl Iterator<Item = (usize, usize)> {
            (0..m).cartesian_product(0..m).filter(|(i, j)| i != j)
        }

        'minimize: loop {
            break {
                for (i, j) in pairs_iter(m) {
                    if Self::leading_term(&generators[j])
                        .unwrap()
                        .divides(&Self::leading_term(&generators[i]).unwrap())
                    {
                        // `LT(g_i)` is superfluous
                        generators.remove(i);
                        m -= 1;

                        continue 'minimize;
                    }
                }
            };
        }

        'reduce: loop {
            break {
                for (i, j) in pairs_iter(m) {
                    for t in generators[i].terms() {
                        // if a nonleading term of `g_i` is divisible by `LT(g_j)`
                        if Self::leading_term(&generators[j]).unwrap().divides(t) {
                            // we will replace `g` with a reduced equivalent
                            let g = generators.remove(i);

                            // `r` is remainder of `g mod (G \ {g})`
                            let (g_rem, _) = Self::divide(&g, &generators);

                            // replaces `g` with
                            generators.push(g_rem);

                            continue 'reduce;
                        }
                    }
                }
            };
        }
    }
}

#[cfg(test)]
mod tests {
    use std::ops::Add;

    use super::*;
    use crate::{
        core::frac::Frac,
        poly::{elts::Polynomial, elts::Term},
    };

    use super::super::ord::Lex;

    type Poly = Polynomial<f64>;
    type Comp = Computer<f64, Lex>;

    #[allow(unused)]
    macro_rules! pp {
        ($poly:expr) => {
            println!("{} = {}", stringify!($poly), &$poly);
        };
    }

    #[allow(unused)]
    macro_rules! pps {
        ($polys:expr) => {
            for i in 0..($polys.len()) {
                println!("{}[{}] = {}", stringify!($polys), i, &$polys[i]);
            }
        };
    }

    #[allow(unused)]
    #[test]
    fn dbg_stuff() {
        let c = |coef: f64| Poly::from(coef);
        fn_vars! { f64: x y z }

        /* let f = y(1);

        let divs = [
            Poly::zero(),
            y(1) + c(1.0),
        ];

        println!("{} / [{}, {}]", &f, &divs[0], &divs[1]);

        let (r, q) = Computer::<Lex, f64>::divide(&f, &divs);

        pp!(r);
        pp!(q[0]);
        pp!(q[1]);

        pp!(&q[0] * &divs[0]);
        pp!(&q[1] * &divs[1]); */

        let f = x(2) + x(3) + x(1) + x(2) + x(0) + x(3);

        pp!(f);
    }

    #[test]
    fn division() {
        let poly_iter = polys(2, 3, 1)
            .into_iter()
            .cartesian_product(polys(2, 1, 2).into_iter().cartesian_product(polys(2, 1, 2)));

        for (f, (g1, g2)) in poly_iter {
            // println!("{} / [{}, {}]", &f, &g1, &g2);
            test_result_equality(&f, &[g1, g2]);
        }
    }

    fn polys(vars: i32, max_deg: u8, max_coef: i32) -> Vec<Polynomial<f64>> {
        (0..vars)
            .map(|_| (0..=max_deg).rev())
            .multi_cartesian_product()
            .map(MultiDegree::from_vec)
            .filter(|mdeg| mdeg.total_deg() <= max_deg)
            .map(|mdeg| {
                (0..=max_coef)
                    .map(f64::from)
                    .map(move |coef| Term::new(coef, mdeg.clone()))
            })
            .multi_cartesian_product()
            .map(|v| Poly::new_unchecked(v.into_iter().filter(|t| !t.is_zero()).collect()))
            .collect()
    }

    #[cfg(test)]
    fn poly_assert_eq(f1: &Poly, f2: &Poly) {
        assert_eq!(Comp::sort_terms(f1).terms, Comp::sort_terms(f2).terms);
    }

    #[cfg(test)]
    fn test_result_equality(f: &Poly, g: &[Poly]) {
        let (r, q) = Comp::divide(f, g);
        let f2 = q.iter().zip_eq(g).map(|(qi, gi)| qi * gi).fold(r, Add::add);
        poly_assert_eq(f, &f2);
    }

    #[test]
    fn reduction() {
        let q = |a, b| Polynomial::from(Frac::new(a, b));
        let c = |coef| q(coef, 1);
        fn_vars! { Frac: x y }
        type CompQ = Computer<Frac, Lex>;

        let f = c(5) * x(4) * y(3) + c(2) * x(2) * y(1) + c(3) * x(1) * y(2) + y(2) + c(3);
        let g = c(8) * x(5) * y(2) + x(3) + c(3) * x(2) * y(2) + y(4) + c(6);

        pp!(f);
        pp!(g);

        println!();

        if let Some(h) = CompQ::try_reduce(&f, &g) {
            pp!(h);
        } else {
            println!("could not reduce");
        }
    }

    #[test]
    fn buchberger() {
        let q = |a, b| Polynomial::from(Frac::new(a, b));
        let c = |coef| q(coef, 1);
        fn_vars! { Frac: x y }
        type CompQ = Computer<Frac, Lex>;

        fn print_buchberger(g: &[Polynomial<Frac>]) {
            let mut g = Vec::from(g);

            println!("initial:");
            pps!(g);

            CompQ::buchberger_algorithm(&mut g);

            println!("\nresult:");
            pps!(g);
            println!("--------------------");
        }

        print_buchberger(&[x(3) * y(1) - x(1) * y(2) + c(1), x(2) * y(2) - y(3) - c(1)]);

        let mut g = vec![x(3) * y(1) - x(1) * y(2) + c(1), x(2) * y(2) - y(3) - c(1)];

        println!("initial:");
        pps!(g);

        CompQ::buchberger_algorithm(&mut g);

        println!("\nresult:");
        pps!(g);

        println!("--------------------");

        let mut g = vec![
            x(2) + x(1) * y(5) + y(4),
            x(1) * y(6) - x(1) * y(3) + y(5) - y(2),
            x(1) * y(5) - x(1) * y(2),
        ];

        println!("initial:");
        pps!(g);

        CompQ::buchberger_algorithm(&mut g);

        println!("\nresult:");
        pps!(g);
    }
}
