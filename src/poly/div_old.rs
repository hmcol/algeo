use std::marker::PhantomData;

use itertools::Itertools;

use crate::core::num::Field;

use super::ord::MonomialOrder;
use super::{MDeg, Polynomial, Term};

/// Handler for computations related to finding Gröbner bases.
///
/// The field `F` is the field of the polynomials it works over and the monomial order `O` is what it uses to determine leading terms
pub struct Computer<F: Field, O: MonomialOrder> {
    _marker: PhantomData<(F, O)>,
}

impl<F: Field, O: MonomialOrder> Computer<F, O> {
    /// Returns the polynomial with the same terms as `f`, sorted by the monomial order `Self::O`.
    pub fn sort_terms(f: &Polynomial<F>) -> Polynomial<F> {
        Polynomial::new_unchecked(
            f.terms()
                .sorted_by(|s, t| O::cmp(&s.mdeg, &t.mdeg))
                .cloned()
                .collect(),
        )
    }

    /// Returns the leading term of `f`, under the monomial order `Self::O`.
    /// 
    /// Returns `None` if `f` has no terms.
    /// 
    /// Ideally, should also return none if it only has zero terms; but also we would ideally never have zero terms, so go figure.
    pub fn leading_term(f: &'_ Polynomial<F>) -> Option<&'_ Term<F>> {
        f.terms().max_by(|s, t| O::cmp(&s.mdeg, &t.mdeg))
    }

    /// Performs general polynomial division on `f` with the divisors `divs`
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
            // maybe not necessary anymore? not sure
            // hack to ignore zero coefficients
            // - will not work in case of floating-point error
            // - should be feature of Polys to cull zeros
            /* if lt_f.coef == F::ZERO {
                f.terms.pop();
                continue;
            } */

            // f still has (nonzero) terms

            for (g, q) in divs.iter().zip(quotients.iter_mut()) {
                if let Some(lt_g) = Self::leading_term(&g) {
                    if lt_g.divides(&lt_f) {
                        // case 1: LT(f) is divisible by some LT(g_i)
                        // - use LT(g_i) to reduce the degree of f

                        // a_i is the term such that LT(f) = a_i * LT(g_i)
                        let a = &(&lt_f / lt_g);
                        // add a_i to g_i's quotient, q_i
                        *q += a;
                        // eliminate LT(f) with a_i * LT(g_i)
                        f -= a * g;

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

                remainder += lt_f;
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
        if s.coef == F::ZERO && t.coef == F::ZERO {
            return Term::zero();
        }

        Term::monic(MDeg(
            s.mdeg
                .degs()
                .zip(t.mdeg.degs())
                .map(|(a, b)| a.max(b))
                .cloned()
                .collect(),
        ))
    }

    /// Returns a combination of `f` and `g` with their leading terms cancelled.
    /// 
    /// In particular, if $M$ is the monic lcm of $\operatorname{LT}(f)$ and $\operatorname{LT}(g)$, `reduce(f, g)` returns
    /// 
    /// $$ S(f, g) = \frac{M}{\operatorname{LT}(f)} f - \frac{M}{\operatorname{LT}(g)} g. $$ 
    pub fn reduce(f: &Polynomial<F>, g: &Polynomial<F>) -> Option<Polynomial<F>> {
        let lt_f = Self::leading_term(f)?;
        let lt_g = Self::leading_term(g)?;
        let m = Self::monic_lcm(&lt_f, &lt_g);
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use std::ops::Add;

    use super::*;
    use crate::poly::{MDeg, Polynomial, Term};

    use super::super::ord::Lex;

    type Poly = Polynomial<f64>;

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
        use crate::poly::{x, y, z};

        let c = |coef| Poly::constant(coef);
        let x = |d| Poly::from(x::<f64>(d));
        let y = |d| Poly::from(y::<f64>(d));
        let z = |d| Poly::from(z::<f64>(d));

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
        const VARS: i32 = 2;
        const MAX_DEG: i8 = 2;
        const MAX_COEF: i32 = 1;

        let poly_iter = (0..VARS)
            .map(|_| (0..=MAX_DEG).rev())
            .multi_cartesian_product()
            .map(MDeg::from_vec)
            .filter(|mdeg| mdeg.total_deg() <= MAX_DEG)
            .map(|mdeg| {
                (0..=MAX_COEF)
                    .map(f64::from)
                    .map(move |coef| Term::new(coef, mdeg.clone()))
            })
            .multi_cartesian_product()
            .map(|v| Poly::new_unchecked(v.into_iter().filter(|t| t.coef != 0.0).collect()));

        for (f, (g1, g2)) in poly_iter
            .clone()
            .cartesian_product(poly_iter.clone().cartesian_product(poly_iter))
        {
            // println!("{} / [{}, {}]", &f, &g1, &g2);
            test_result_equality(&f, &[g1, g2]);
        }
    }

    #[cfg(test)]
    fn poly_assert_eq(f1: &Poly, f2: &Poly) {
        assert_eq!(
            Computer::<f64, Lex>::sort_terms(f1).terms,
            Computer::<f64, Lex>::sort_terms(f2).terms
        );
    }

    #[cfg(test)]
    fn test_result_equality(f: &Poly, g: &[Poly]) {
        let (r, q) = Computer::<f64, Lex>::divide(f, g);

        let f2 = q.iter().zip_eq(g).map(|(qi, gi)| qi * gi).fold(r, Add::add);

        poly_assert_eq(f, &f2);
    }
}
