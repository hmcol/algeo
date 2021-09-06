use itertools::Itertools;

use crate::core::num::Field;

use super::ord::MonomialOrder;
use super::{Polynomial, Term};

pub fn sort_terms<MO: MonomialOrder, F: Field>(f: &Polynomial<F>) -> Polynomial<F> {
    Polynomial::from_vec(
        f.terms()
            .sorted_by(|s, t| MO::cmp(&s.mdeg, &t.mdeg))
            .cloned()
            .collect(),
    )
}

pub fn leading_term<'a, MO: MonomialOrder, F: Field>(f: &'a Polynomial<F>) -> Option<&'a Term<F>> {
    f.terms().max_by(|s, t| MO::cmp(&s.mdeg, &t.mdeg))
}


fn divide<MO: MonomialOrder, F: Field>(f: &Polynomial<F>, divs: &[Polynomial<F>]) -> (Polynomial<F>, Vec<Polynomial<F>>) {
    let m = divs.len();
        let mut quotients = vec![Polynomial::<F>::zero(); m];
        let mut remainder = Polynomial::<F>::zero();

        let mut f = f.clone();

        'outer: loop {
            if let Some(lt_f) = leading_term::<MO, F>(&f).cloned() {
                // hack to ignore zero coefficients
                // - will not work in case of floating-point error
                // - should be feature of polynomials to cull zeros
                if lt_f.coef == F::ZERO {
                    f.terms.pop();
                    continue;
                }

                // f still has (nonzero) terms

                for (g, q) in divs.iter().zip(quotients.iter_mut()) {
                    if let Some(lt_g) = leading_term::<MO, F>(&g) {
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
                f.terms.sort_by(|s, t| MO::cmp(&s.mdeg, &t.mdeg));

                // ofter sorting, last term of f is the leading term, LT(f)
                if let Some(lt_f) = f.terms.pop() {
                    // first

                    remainder += lt_f;
                }
            } else {
                // f has no more (nonzero) terms; division is done
                break;
            }
        }

        // return
        (remainder, quotients)
}


#[cfg(test)]
mod tests {
    use std::ops::Add;

    use super::*;
    use crate::poly::*;

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
    fn dbg_stuff() {
        let c = |coef| Poly::constant(coef);
        let x = |d| Poly::from(x::<f64>(d));
        let y = |d| Poly::from(y::<f64>(d));
        let z = |d| Poly::from(z::<f64>(d));


        let f = <Poly as Add<Poly>>::add(x(1), y(1));
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
            .map(|v| Poly::from_vec(v.into_iter().filter(|t| t.coef != 0.0).collect()));


        for (f, (g1, g2)) in poly_iter
            .clone()
            .cartesian_product(poly_iter.clone().cartesian_product(poly_iter))
        {
            // println!("{} / [{}, {}]", &f, &g1, &g2);
            test_result_equality::<Lex, f64>(&f, &[g1, g2]);
        }
    }

    #[cfg(test)]
    fn poly_assert_eq<MO: MonomialOrder, F: Field>(f1: &Polynomial<F>, f2: &Polynomial<F>) {
        assert_eq!(sort_terms::<MO, F>(f1).terms, sort_terms::<MO, F>(f2).terms);
    }

    #[cfg(test)]
    fn test_result_equality<MO: MonomialOrder, F: Field>(
        f: &Polynomial<F>,
        g: &[Polynomial<F>],
    ) {
        let (r, q) = divide::<MO, F>(f, g);

        let f2 = q.iter().zip_eq(g).map(|(qi, gi)| qi * gi).fold(r, Add::add);

        poly_assert_eq::<MO, F>(f, &f2);
    }
}