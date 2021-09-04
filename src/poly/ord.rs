use std::{cmp::Ordering, marker::PhantomData};

use itertools::Itertools;

use super::{MDeg, Polynomial as Poly, Term};
use crate::core::num::Field;

/// The [Lexicographic Order](https://w.wiki/3zwi) on multidegrees.
///
/// if a != b, compares first unequal degrees from the left
///
/// e.g., a < b iff ∃k s.t. a_k < b_k and a_i = b_i, for i = 0,...,k-1
pub fn lex(a: &MDeg, b: &MDeg) -> Ordering {
    for (deg_a, deg_b) in a.degs().zip(b.degs()) {
        match deg_a.cmp(deg_b) {
            Ordering::Equal => continue,
            lt_or_gt => return lt_or_gt,
        }
    }
    grad(a, b)
}

/// The Reverse [Lexicographic Order](https://w.wiki/3zwi) order on
/// multidegrees.
///
/// This runs the lexicographic order with the indices reversed; not to be
/// confused with simply calling `Ordering::reverse` on the result of [`lex`].
pub fn revlex(a: &MDeg, b: &MDeg) -> Ordering {
    match a.len().cmp(&b.len()) {
        Ordering::Equal => {
            for (deg_a, deg_b) in a.degs().zip(b.degs()).rev() {
                match deg_a.cmp(deg_b) {
                    Ordering::Equal => continue,
                    lt_or_gt => return lt_or_gt,
                }
            }
        }
        lt_or_gt => return lt_or_gt,
    }
    Ordering::Equal
}

/// The graded order on multidegrees.
///
/// Simply compares the total degrees.
///
/// This is the usual grading on a univariate polynomial ring.
///
/// Important note: this is not a 'monomial order' as it is not antisymmetric;
/// should probably be moved or hidden to avoid confusion
pub fn grad(a: &MDeg, b: &MDeg) -> Ordering {
    a.total_deg().cmp(&b.total_deg())
}

/// The [Graded Lexicographic Order](https://w.wiki/3zwp) on multidegrees.
///
/// applies the graded order; if equal, applies lexicographic
pub fn grlex(a: &MDeg, b: &MDeg) -> Ordering {
    match grad(a, b) {
        Ordering::Equal => lex(a, b),
        lt_or_gt => lt_or_gt,
    }
}

/// The [Graded Reverse Lexicographic Order](https://w.wiki/3zwq) on
/// multidegrees.
///
/// applies the graded order; if equal, applies reverse lexicographic with the result negated
pub fn grevlex(a: &MDeg, b: &MDeg) -> Ordering {
    match grad(a, b) {
        Ordering::Equal => revlex(a, b).reverse(),
        lt_or_gt => lt_or_gt,
    }
}

/// More possible orders:
/// - https://en.wikipedia.org/wiki/Monomial_order
/// - https://faculty.math.illinois.edu/Macaulay2/doc/Macaulay2-1.15/share/doc/Macaulay2/Macaulay2Doc/html/_monomial_sporderings.html
///
struct _PlaceHolder;

/* enum Ordered<T> {
    Lex(T),
    RevLex(T),
    Grad(T),
    GrLex(T),
    GrevLex(T),
}

impl<T> Ordered<T> {
    pub fn inner<'a>(&'a self) -> &'a T {
        match self {
            Ordered::Lex(x) => x,
            Ordered::RevLex(x) => x,
            Ordered::Grad(x) => x,
            Ordered::GrLex(x) => x,
            Ordered::GrevLex(x) => x,
        }
    }
}

impl PartialEq for Ordered<MDeg> {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::Lex(a), Self::Lex(b)) => a == b,
            (Self::RevLex(a), Self::RevLex(b)) => a == b,
            (Self::Grad(a), Self::Grad(b)) => a == b,
            (Self::GrLex(a), Self::GrLex(b)) => a == b,
            (Self::GrevLex(a), Self::GrevLex(b)) => a == b,
            _ => false,
        }
    }
}

impl Eq for Ordered<MDeg> {}

impl PartialOrd for Ordered<MDeg> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        match (self, other) {
            (Self::Lex(a), Self::Lex(b)) => Some(lex(a, b)),
            (Self::RevLex(a), Self::RevLex(b)) => Some(revlex(a, b)),
            (Self::Grad(a), Self::Grad(b)) => Some(grad(a, b)),
            (Self::GrLex(a), Self::GrLex(b)) => Some(grlex(a, b)),
            (Self::GrevLex(a), Self::GrevLex(b)) => Some(grevlex(a, b)),
            _ => None,
        }
    }
}

impl Ord for Ordered<MDeg> {
    fn cmp(&self, other: &Self) -> Ordering {
        // should handle monomial ordering better; this is stinky
        self.partial_cmp(other).expect(
            "illegal comparison of ordered multidegrees of differing \
            variants; this error implies a flaw in the algeo crate",
        )
    }
}

impl<F: Field> Ordered<Poly<F>> {
    pub fn leading_term<'a>(&'a self) -> Ordered<&'a Term<F>> {
        todo!()
    }
} */

#[derive(Clone, Copy)]
struct DivisionComputer<F: Field> {
    order: fn(&MDeg, &MDeg) -> Ordering,
    _field_marker: PhantomData<F>,
}

impl<F: Field> DivisionComputer<F> {
    pub fn new(order: fn(&MDeg, &MDeg) -> Ordering) -> Self {
        DivisionComputer {
            order: order,
            _field_marker: PhantomData,
        }
    }

    pub fn divide(&self, f: &Poly<F>, divs: &[Poly<F>]) -> (Poly<F>, Vec<Poly<F>>) {
        let m = divs.len();
        let mut quotients = vec![Poly::<F>::zero(); m];
        let mut remainder = Poly::<F>::zero();

        let mut f = f.clone();

        'outer: loop {
            if let Some(lt_f) = self.leading_term(&f).cloned() {
                // hack to ignore zero coefficients
                // - will not work in case of floating-point error
                // - should be feature of polynomials to cull zeros
                if lt_f.coef == F::ZERO {
                    f.terms.pop();
                    continue;
                }

                // f still has (nonzero) terms

                for (g, q) in divs.iter().zip(quotients.iter_mut()) {
                    if let Some(lt_g) = self.leading_term(&g) {
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
                f.terms.sort_by(|s, t| (self.order)(&s.mdeg, &t.mdeg));

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

    pub fn sort_terms(&self, f: &Poly<F>) -> Poly<F> {
        Poly::from_vec(
            f.terms()
                .sorted_by(|s, t| (self.order)(&s.mdeg, &t.mdeg))
                .cloned()
                .collect(),
        )
    }

    pub fn leading_term<'a>(&self, f: &'a Poly<F>) -> Option<&'a Term<F>> {
        f.terms().max_by(|s, t| (self.order)(&s.mdeg, &t.mdeg))
    }
}

#[cfg(test)]
mod tests {
    use itertools::Itertools;
    use std::cmp::Ordering;

    use super::DivisionComputer;
    use crate::poly::*;

    type Poly = Polynomial<f64>;

    macro_rules! mdeg {
        ($( $deg:expr ),* $(,)?) => {
            &MDeg::from_vec(vec![ $( $deg ),* ])
        };
    }

    fn dbg_suite(ord: fn(&MDeg, &MDeg) -> Ordering) {
        let d = |s: &[i8]| format!("{:?}{:?}{:?}", s[0], s[1], s[2]);

        let c = |o: Ordering| match o {
            Ordering::Less => '<',
            Ordering::Equal => '=',
            Ordering::Greater => '>',
        };

        let iter = (0..3).map(|_| 0..3).multi_cartesian_product();

        for (ref a, ref b) in iter.clone().cartesian_product(iter) {
            let mdeg_a = &MDeg::from_vec(a.clone());
            let mdeg_b = &MDeg::from_vec(b.clone());

            if super::grad(mdeg_a, mdeg_b).is_eq() {
                println!("{} {} {}", d(a), c(ord(mdeg_a, mdeg_b)), d(b));
            }
        }
    }

    #[test]
    fn test_lex() {
        use super::lex;

        assert_eq!(lex(mdeg![0, 0, 0], mdeg![0, 0, 0]), Ordering::Equal);
        assert_eq!(lex(mdeg![1, 0, 0], mdeg![1, 0, 0]), Ordering::Equal);
        assert_eq!(lex(mdeg![1, 2, 3], mdeg![1, 2, 3]), Ordering::Equal);
        assert_eq!(lex(mdeg![0, 0, 1], mdeg![0, 0, 1]), Ordering::Equal);

        assert_eq!(lex(mdeg![0, 1, 0], mdeg![1, 0, 1]), Ordering::Less);
        assert_eq!(lex(mdeg![0, 0, 1], mdeg![1, 0, 0]), Ordering::Less);

        assert_eq!(lex(mdeg![1, 0, 0], mdeg![0, 1, 0]), Ordering::Greater);
        assert_eq!(lex(mdeg![2, 2, 0], mdeg![0, 0, 1]), Ordering::Greater);

        // dbg_suite(lex);
    }

    #[test]
    fn test_revlex() {
        use super::revlex;

        assert_eq!(revlex(mdeg![0, 0, 0], mdeg![0, 0, 0]), Ordering::Equal);
        assert_eq!(revlex(mdeg![1, 0, 0], mdeg![1, 0, 0]), Ordering::Equal);
        assert_eq!(revlex(mdeg![1, 2, 3], mdeg![1, 2, 3]), Ordering::Equal);
        assert_eq!(revlex(mdeg![0, 0, 1], mdeg![0, 0, 1]), Ordering::Equal);

        assert_eq!(revlex(mdeg![1, 0, 0], mdeg![0, 1, 0]), Ordering::Less);
        assert_eq!(revlex(mdeg![0, 1, 0], mdeg![1, 0, 1]), Ordering::Less);
        assert_eq!(revlex(mdeg![2, 2, 0], mdeg![0, 0, 1]), Ordering::Less);

        assert_eq!(revlex(mdeg![0, 0, 1], mdeg![1, 0, 0]), Ordering::Greater);

        // suite(revlex);
    }

    #[test]
    fn test_grad() {
        use super::grad;

        assert_eq!(grad(mdeg![0, 0, 0], mdeg![0, 0, 0]), Ordering::Equal);
        assert_eq!(grad(mdeg![1, 0, 0], mdeg![1, 0, 0]), Ordering::Equal);
        assert_eq!(grad(mdeg![1, 2, 3], mdeg![1, 2, 3]), Ordering::Equal);
        assert_eq!(grad(mdeg![0, 0, 1], mdeg![0, 0, 1]), Ordering::Equal);
        assert_eq!(grad(mdeg![1, 0, 0], mdeg![0, 1, 0]), Ordering::Equal);

        assert_eq!(grad(mdeg![2, 2, 0], mdeg![1, 1, 5]), Ordering::Less);
        assert_eq!(grad(mdeg![1, 0, 0], mdeg![0, 1, 1]), Ordering::Less);

        assert_eq!(grad(mdeg![2, 2, 0], mdeg![0, 0, 1]), Ordering::Greater);
        assert_eq!(grad(mdeg![0, 0, 2], mdeg![1, 0, 0]), Ordering::Greater);
        assert_eq!(grad(mdeg![0, 3, 0], mdeg![1, 0, 1]), Ordering::Greater);

        // dbg_suite(grad);
    }

    #[test]
    fn test_grlex() {
        use super::grlex;

        assert_eq!(grlex(mdeg![0, 0, 0], mdeg![0, 0, 0]), Ordering::Equal);
        assert_eq!(grlex(mdeg![1, 0, 0], mdeg![1, 0, 0]), Ordering::Equal);
        assert_eq!(grlex(mdeg![1, 2, 3], mdeg![1, 2, 3]), Ordering::Equal);
        assert_eq!(grlex(mdeg![0, 0, 1], mdeg![0, 0, 1]), Ordering::Equal);

        assert_eq!(grlex(mdeg![2, 2, 0], mdeg![1, 1, 5]), Ordering::Less);
        assert_eq!(grlex(mdeg![0, 1, 1], mdeg![1, 1, 0]), Ordering::Less);
        assert_eq!(grlex(mdeg![0, 3, 0], mdeg![1, 0, 2]), Ordering::Less);

        assert_eq!(grlex(mdeg![1, 0, 0], mdeg![0, 1, 0]), Ordering::Greater);
        assert_eq!(grlex(mdeg![2, 2, 0], mdeg![0, 0, 1]), Ordering::Greater);
        assert_eq!(grlex(mdeg![0, 0, 2], mdeg![1, 0, 0]), Ordering::Greater);
        assert_eq!(grlex(mdeg![3, 0, 0], mdeg![1, 0, 2]), Ordering::Greater);

        // dbg_suite(grlex);
    }

    #[test]
    fn test_grevlex() {
        use super::grevlex;

        assert_eq!(grevlex(mdeg![0, 0, 0], mdeg![0, 0, 0]), Ordering::Equal);
        assert_eq!(grevlex(mdeg![1, 0, 0], mdeg![1, 0, 0]), Ordering::Equal);
        assert_eq!(grevlex(mdeg![1, 2, 3], mdeg![1, 2, 3]), Ordering::Equal);
        assert_eq!(grevlex(mdeg![0, 0, 1], mdeg![0, 0, 1]), Ordering::Equal);
        

        assert_eq!(grevlex(mdeg![2, 2, 0], mdeg![1, 1, 5]), Ordering::Less);
        assert_eq!(grevlex(mdeg![1, 0, 0], mdeg![0, 1, 1]), Ordering::Less);

        assert_eq!(grevlex(mdeg![2, 2, 0], mdeg![0, 0, 1]), Ordering::Greater);
        assert_eq!(grevlex(mdeg![0, 0, 2], mdeg![1, 0, 0]), Ordering::Greater);
        assert_eq!(grevlex(mdeg![0, 3, 0], mdeg![1, 0, 1]), Ordering::Greater);
        assert_eq!(grevlex(mdeg![1, 0, 0], mdeg![0, 1, 0]), Ordering::Greater);

        // dbg_suite(grevlex);
    }

    // #[test]
    fn dbg_grlex_vs_grevlex() {
        let mut vecs = Vec::new();

        const D: i8 = 5;

        for x in 0..D {
            for y in 0..(D - x) {
                for z in 0..(D - x - y) {
                    vecs.push(MDeg::from_vec(vec![x, y, z]));
                }
            }
        }

        for (a, b) in vecs.iter().cartesian_product(&vecs) {
            let grlex = super::grlex(a, b);
            let grevlex = super::grevlex(a, b);

            if grlex != grevlex {
                println!(
                    "cmp {} {} => grlex: {:<10} grevlex: {:<10}",
                    a,
                    b,
                    format!("{:?}", grlex),
                    format!("{:?}", grevlex)
                );
            }
        }
    }

    fn print_poly(head: &str, f: &Poly) {
        println!("{}", format!("{} {}", head, f));
    }

    macro_rules! pp {
        ($poly:expr) => {
            println!("{} = {}", stringify!($poly), &$poly);
        };
    }

    macro_rules! pps {
        ($polys:expr) => {
            for i in 0..($polys.len()) {
                println!("{}[{}] = {}", stringify!($polys), i, &$polys[i]);
            }
        };
    }

    #[test]
    fn division() {
        let c = |coef| Poly::from(Term::<f64>::from(coef));
        let x = |d| Poly::from(x::<f64>(d));
        let y = |d| Poly::from(y::<f64>(d));
        // let z = |d| Poly::from(z::<f64>(d));
        
        let dc = DivisionComputer::<f64>::new(lex);
        let div = |f, divs| dc.divide(f, divs);

        let f = c(1.0) * x(3) * y(3) + c(3.0) * x(2) * y(4);
        let g = [c(1.0) * x(1) * y(4)];

        pp!(f);
        pps!(g);
        println!();

        let (r, q) = div(&f, &g);
        pp!(r);
        pps!(q);
        println!();

        let f2 = &q[0] * &g[0] + r;
        pp!(f2);

        println!("\n------------------------------\n");

        let f = x(2) + x(1) - y(2) + y(1);
        let g = [x(1) * y(1) + c(1.0), x(1) + y(1)];

        pp!(f);
        pps!(g);
        println!();

        let (r, q) = div(&f, &g);
        pp!(r);
        pps!(q);
        println!();

        let f2 = &q[0] * &g[0] + &q[1] * &g[1] + r;
        pp!(f2);

        println!("\n------------------------------\n");

        let f = x(2) + x(1) - y(2) + y(1);
        let g = [x(1) + y(1), x(1) * y(1) + c(1.0)];

        pp!(f);
        pps!(g);
        println!();

        let (r, q) = div(&f, &g);
        pp!(r);
        pps!(q);
        println!();

        let f2 = &q[0] * &g[0] + &q[1] * &g[1] + r;
        pp!(f2);
    }
}
