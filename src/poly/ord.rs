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
///
/// Currently gives the wrong answer if both are equal up to some index, and
/// one has only remaining zeros (i.e., a lack of entries with greater index).
/// In this case, the multidegree with only zeros remaining should be less.
pub fn lex(ref a: &MDeg, b: &MDeg) -> Ordering {
    for ((idx_a, deg_a), (idx_b, deg_b)) in a.degs().zip(b.degs()) {
        match idx_a.cmp(idx_b) {
            // Both `a` and `b` have nonzero entries at the given index.
            Ordering::Equal => match deg_a.cmp(deg_b) {
                Ordering::Equal => continue,
                lt_or_gt => return lt_or_gt,
            },
            // Either `a` or `b` has a zero at an index where the other has a
            // nonzero degree. The index `idx_a` (resp. `idx_b`) is the first
            // nonzero index of `a` (resp. `b`) where the two differ. The
            // greater of the two has a zero where the other has a positive
            // degree. Therefore, the lesser multidegree is the one with the
            // greater index, here.
            lt_or_gt => return lt_or_gt.reverse(),
        }
    }
    grad(a, b)
}

/// The Reverse [Lexicographic Order](https://w.wiki/3zwi) order on
/// multidegrees.
///
/// This runs the lexicographic order with the indices reversed; not to be
/// confused with simply calling `Ordering::reverse` on the result of [`lex`].
///
/// Same issue as lexicographic.
///
/// Likely some way to avoid this repetition or a better approach in general.
pub fn revlex(a: &MDeg, b: &MDeg) -> Ordering {
    for ((idx_a, deg_a), (idx_b, deg_b)) in a.degs().zip(b.degs()).rev() {
        match idx_a.cmp(idx_b) {
            Ordering::Equal => match deg_a.cmp(deg_b) {
                Ordering::Equal => continue,
                lt_or_gt => return lt_or_gt,
            },
            lt_or_gt => return lt_or_gt.reverse(),
        }
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
struct PlaceHolder;

#[cfg(test)]
mod tests {
    use std::cmp::Ordering;

    use itertools::Itertools;

    use crate::poly::*;

    macro_rules! mdeg {
        ($( $deg:expr ),* $(,)?) => {
            &MDeg::from_slice(&[ $( $deg ),* ])
        };
    }

    #[test]
    fn lex() {
        use super::lex;

        dbg!(mdeg![3, 0 , 2, 0, 0, 1, 4, 0, 0, 9]);

        let d = |s: &[i8]| format!("{:?}{:?}{:?}", s[0], s[1], s[2]);

        let print = |(ref a, ref b): (Vec<i8>, Vec<i8>)| {
            let mdeg_a = &MDeg::from_slice(a);
            let mdeg_b = &MDeg::from_slice(b);

            println!("{} ~ {} => {:?}", d(a), d(b), lex(mdeg_a, mdeg_b));
        };

        assert_eq!(lex(mdeg![1], mdeg![1]), Ordering::Equal);
        assert_eq!(lex(mdeg![1, 2, 3], mdeg![1, 2, 3]), Ordering::Equal);

        assert_eq!(lex(mdeg![1], mdeg![0, 1]), Ordering::Greater);

        assert_eq!(lex(mdeg![0, 3, 2], mdeg![0, 3, 3]), Ordering::Less);

        let iter = (0..=2).map(|_| 0..=1i8).multi_cartesian_product();

        iter.clone().cartesian_product(iter).for_each(print);
    }
}
