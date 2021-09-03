use std::cmp::Ordering;

use super::{MDeg, Polynomial, Term};
use crate::core::num::Field;

/// enum for various possible monomial orders on multidegrees
pub enum MonomialOrder {
    Lex,
    RevLex,
    Grad,
    GrLex,
    GRevLex,
}

/// applies the given monomial order to a and b
pub fn cmp(order: MonomialOrder, a: &MDeg, b: &MDeg) -> Ordering {
    match order {
        MonomialOrder::Lex => lex(a, b),
        MonomialOrder::RevLex => revlex(a, b),
        MonomialOrder::Grad => grad(a, b),
        MonomialOrder::GrLex => grlex(a, b),
        MonomialOrder::GRevLex => grevlex(a, b),
    }
}

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




#[cfg(test)]
mod tests {
    use itertools::Itertools;
    use std::cmp::Ordering;

    use crate::poly::*;

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

        dbg_suite(grlex);
    }

    #[test]
    fn test_grevlex() {
        use super::grevlex;

        assert_eq!(grevlex(mdeg![0, 0, 0], mdeg![0, 0, 0]), Ordering::Equal);
        assert_eq!(grevlex(mdeg![1, 0, 0], mdeg![1, 0, 0]), Ordering::Equal);
        assert_eq!(grevlex(mdeg![1, 2, 3], mdeg![1, 2, 3]), Ordering::Equal);
        assert_eq!(grevlex(mdeg![0, 0, 1], mdeg![0, 0, 1]), Ordering::Equal);
        assert_eq!(grevlex(mdeg![1, 0, 0], mdeg![0, 1, 0]), Ordering::Equal);

        assert_eq!(grevlex(mdeg![2, 2, 0], mdeg![1, 1, 5]), Ordering::Less);
        assert_eq!(grevlex(mdeg![1, 0, 0], mdeg![0, 1, 1]), Ordering::Less);

        assert_eq!(grevlex(mdeg![2, 2, 0], mdeg![0, 0, 1]), Ordering::Greater);
        assert_eq!(grevlex(mdeg![0, 0, 2], mdeg![1, 0, 0]), Ordering::Greater);
        assert_eq!(grevlex(mdeg![0, 3, 0], mdeg![1, 0, 1]), Ordering::Greater);

        // dbg_suite(grevlex);
    }

    #[test]
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
}
