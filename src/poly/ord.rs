use std::cmp::Ordering;

use super::MDeg;

pub trait MonomialOrder {
    fn cmp(a: &MDeg, b: &MDeg) -> Ordering;
}

/// The [Lexicographic Order](https://w.wiki/3zwi) on multidegrees.
///
/// if a != b, compares first unequal degrees from the left
///
/// e.g., a < b iff ∃k s.t. a_k < b_k and a_i = b_i, for i = 0,...,k-1
pub struct Lex;

impl MonomialOrder for Lex {
    fn cmp(a: &MDeg, b: &MDeg) -> Ordering {
        for (deg_a, deg_b) in a.degs().zip(b.degs()) {
            match deg_a.cmp(deg_b) {
                Ordering::Equal => continue,
                lt_or_gt => return lt_or_gt,
            }
        }
        grad(a, b)
    }
}

/// The Reverse [Lexicographic Order](https://w.wiki/3zwi) order on
/// multidegrees.
///
/// This runs the lexicographic order with the indices reversed; not to be
/// confused with simply calling `Ordering::reverse` on the result of [`lex`].
pub struct RevLex;

impl MonomialOrder for RevLex {
    fn cmp(a: &MDeg, b: &MDeg) -> Ordering {
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
pub struct GrLex;

impl MonomialOrder for GrLex {
    fn cmp(a: &MDeg, b: &MDeg) -> Ordering {
        match grad(a, b) {
            Ordering::Equal => Lex::cmp(a, b),
            lt_or_gt => lt_or_gt,
        }
    }
}

/// The [Graded Reverse Lexicographic Order](https://w.wiki/3zwq) on
/// multidegrees.
///
/// applies the graded order; if equal, applies reverse lexicographic with the result negated
pub struct GRevLex;

impl MonomialOrder for GRevLex {
    fn cmp(a: &MDeg, b: &MDeg) -> Ordering {
        match grad(a, b) {
            Ordering::Equal => RevLex::cmp(a, b).reverse(),
            lt_or_gt => lt_or_gt,
        }
    }
}

/// More possible orders:
/// - https://en.wikipedia.org/wiki/Monomial_order
/// - https://faculty.math.illinois.edu/Macaulay2/doc/Macaulay2-1.15/share/doc/Macaulay2/Macaulay2Doc/html/_monomial_sporderings.html
///
struct _PlaceHolder;

#[cfg(test)]
mod order_tests {
    use itertools::Itertools;
    use std::cmp::Ordering;

    use crate::poly::{MDeg, MonomialOrder};
    use super::*;

    macro_rules! mdeg {
        ($( $deg:expr ),* $(,)?) => {
            &MDeg::from_vec(vec![ $( $deg ),* ])
        };
    }

    fn _dbg_suite(ord: fn(&MDeg, &MDeg) -> Ordering) {
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

    #[cfg(test)]
    macro_rules! assert_ord_equal {
        ($ord:expr) => {
            assert_eq!($ord, Ordering::Equal)
        };
    }

    #[cfg(test)]
    macro_rules! assert_ord_less {
        ($ord:expr) => {
            assert_eq!($ord, Ordering::Less)
        };
    }

    #[cfg(test)]
    macro_rules! assert_ord_greater {
        ($ord:expr) => {
            assert_eq!($ord, Ordering::Greater)
        };
    }

    #[test]
    fn test_grad() {


        assert_ord_equal!(grad(mdeg![0, 0, 0], mdeg![0, 0, 0]));
        assert_ord_equal!(grad(mdeg![1, 0, 0], mdeg![1, 0, 0]));
        assert_ord_equal!(grad(mdeg![1, 2, 3], mdeg![1, 2, 3]));
        assert_ord_equal!(grad(mdeg![0, 0, 1], mdeg![0, 0, 1]));
        assert_ord_equal!(grad(mdeg![1, 0, 0], mdeg![0, 1, 0]));

        assert_ord_less!(grad(mdeg![2, 2, 0], mdeg![1, 1, 5]));
        assert_ord_less!(grad(mdeg![1, 0, 0], mdeg![0, 1, 1]));

        assert_ord_greater!(grad(mdeg![2, 2, 0], mdeg![0, 0, 1]));
        assert_ord_greater!(grad(mdeg![0, 0, 2], mdeg![1, 0, 0]));
        assert_ord_greater!(grad(mdeg![0, 3, 0], mdeg![1, 0, 1]));

        // dbg_suite(grad);
    }

    #[test]
    fn test_lex() {
        assert_ord_equal!(Lex::cmp(mdeg![0, 0, 0], mdeg![0, 0, 0]));
        assert_ord_equal!(Lex::cmp(mdeg![1, 0, 0], mdeg![1, 0, 0]));
        assert_ord_equal!(Lex::cmp(mdeg![1, 2, 3], mdeg![1, 2, 3]));
        assert_ord_equal!(Lex::cmp(mdeg![0, 0, 1], mdeg![0, 0, 1]));

        assert_ord_less!(Lex::cmp(mdeg![0, 1, 0], mdeg![1, 0, 1]));
        assert_ord_less!(Lex::cmp(mdeg![0, 0, 1], mdeg![1, 0, 0]));

        assert_ord_greater!(Lex::cmp(mdeg![1, 0, 0], mdeg![0, 1, 0]));
        assert_ord_greater!(Lex::cmp(mdeg![2, 2, 0], mdeg![0, 0, 1]));

        // dbg_suite(Lex::cmp);
    }

    #[test]
    fn test_revlex() {
        assert_ord_equal!(RevLex::cmp(mdeg![0, 0, 0], mdeg![0, 0, 0]));
        assert_ord_equal!(RevLex::cmp(mdeg![1, 0, 0], mdeg![1, 0, 0]));
        assert_ord_equal!(RevLex::cmp(mdeg![1, 2, 3], mdeg![1, 2, 3]));
        assert_ord_equal!(RevLex::cmp(mdeg![0, 0, 1], mdeg![0, 0, 1]));

        assert_ord_less!(RevLex::cmp(mdeg![1, 0, 0], mdeg![0, 1, 0]));
        assert_ord_less!(RevLex::cmp(mdeg![0, 1, 0], mdeg![1, 0, 1]));
        assert_ord_less!(RevLex::cmp(mdeg![2, 2, 0], mdeg![0, 0, 1]));

        assert_ord_greater!(RevLex::cmp(mdeg![0, 0, 1], mdeg![1, 0, 0]));

        // suite(RevLex::cmp);
    }

    #[test]
    fn test_grlex() {
        assert_ord_equal!(GrLex::cmp(mdeg![0, 0, 0], mdeg![0, 0, 0]));
        assert_ord_equal!(GrLex::cmp(mdeg![1, 0, 0], mdeg![1, 0, 0]));
        assert_ord_equal!(GrLex::cmp(mdeg![1, 2, 3], mdeg![1, 2, 3]));
        assert_ord_equal!(GrLex::cmp(mdeg![0, 0, 1], mdeg![0, 0, 1]));

        assert_ord_less!(GrLex::cmp(mdeg![2, 2, 0], mdeg![1, 1, 5]));
        assert_ord_less!(GrLex::cmp(mdeg![0, 1, 1], mdeg![1, 1, 0]));
        assert_ord_less!(GrLex::cmp(mdeg![0, 3, 0], mdeg![1, 0, 2]));

        assert_ord_greater!(GrLex::cmp(mdeg![1, 0, 0], mdeg![0, 1, 0]));
        assert_ord_greater!(GrLex::cmp(mdeg![2, 2, 0], mdeg![0, 0, 1]));
        assert_ord_greater!(GrLex::cmp(mdeg![0, 0, 2], mdeg![1, 0, 0]));
        assert_ord_greater!(GrLex::cmp(mdeg![3, 0, 0], mdeg![1, 0, 2]));

        // dbg_suite(GrLex::cmp);
    }

    #[test]
    fn test_grevlex() {
        use super::GRevLex;

        assert_ord_equal!(GRevLex::cmp(mdeg![0, 0, 0], mdeg![0, 0, 0]));
        assert_ord_equal!(GRevLex::cmp(mdeg![1, 0, 0], mdeg![1, 0, 0]));
        assert_ord_equal!(GRevLex::cmp(mdeg![1, 2, 3], mdeg![1, 2, 3]));
        assert_ord_equal!(GRevLex::cmp(mdeg![0, 0, 1], mdeg![0, 0, 1]));

        assert_ord_less!(GRevLex::cmp(mdeg![2, 2, 0], mdeg![1, 1, 5]));
        assert_ord_less!(GRevLex::cmp(mdeg![1, 0, 0], mdeg![0, 1, 1]));

        assert_ord_greater!(GRevLex::cmp(mdeg![2, 2, 0], mdeg![0, 0, 1]));
        assert_ord_greater!(GRevLex::cmp(mdeg![0, 0, 2], mdeg![1, 0, 0]));
        assert_ord_greater!(GRevLex::cmp(mdeg![0, 3, 0], mdeg![1, 0, 1]));
        assert_ord_greater!(GRevLex::cmp(mdeg![1, 0, 0], mdeg![0, 1, 0]));

        // dbg_suite(GRevLex::cmp);
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
            let grlex = GrLex::cmp(a, b);
            let grevlex = GRevLex::cmp(a, b);

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
