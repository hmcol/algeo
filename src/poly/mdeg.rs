use itertools::{EitherOrBoth, Itertools};
use std::cmp::Ordering;
use std::ops::Add;
use xops::binop;

// structs ---------------------------------------------------------------------

/// type for indexing the indeterminates
type I = usize;

/// type of the degrees in a multidegree
///
/// may need to be changed to `u8` as many of the computations are turning out to require nonnegative degrees
type D = u8;

/// Multidegree for a monomial; wraps a `Vec<D>`.
///
/// This is the treatment of multidegrees as (nonnegative) integer tuples.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct MultiDegree(pub Vec<D>);

impl MultiDegree {
    /// Returns an empty multidegree.
    ///
    /// In any mathematical context, you should prefer `MDeg::zero`.
    #[inline]
    pub fn new() -> Self {
        MultiDegree::zero()
    }

    /// Quick fix for sanitizing a multidegree.
    /// 
    /// Should be avoided if possible; prefer clean operations.
    #[inline]
    pub fn trim_zeros(&mut self) {
        if let Some(idx) = self.degs().rposition(|deg| *deg != 0) {
            self.0.truncate(idx + 1);
        }
    }

    /// Returns the multidegree wrapping the given vector, after truncating off
    /// any trailing zeros
    pub fn from_vec(mut vec: Vec<D>) -> Self {
        if let Some(idx) = vec.iter().rposition(|deg| *deg != 0) {
            vec.truncate(idx + 1);
            MultiDegree(vec)
        } else {
            MultiDegree::zero()
        }
    }

    /// returns the empty multidegree
    ///
    /// used for constant terms, i.e., terms without indeterminates
    #[inline]
    pub fn zero() -> Self {
        MultiDegree(Vec::new())
    }

    // Check if `self` is the zero multidegree.
    pub fn is_zero(&self) -> bool {
        self.0.is_empty()
    }

    /// Returns the multidegree with the only nonzero entry of `deg` at `idx`.
    /// 
    /// This is the multidegree of the indeterminate variable `x_idx^deg`.
    pub fn var(idx: I, deg: D) -> Self {
        if deg == 0 {
            MultiDegree::zero()
        } else {
            let mut degs = vec![0; idx];
            degs.push(deg);
            MultiDegree(degs)
        }
    }

    /// returns iterator over the individual degree components
    ///
    /// Should replace return with in-house iterator struct.
    #[inline]
    pub fn degs(&self) -> std::slice::Iter<D> {
        self.0.iter()
    }

    /// returns mutable iterator over the individual degree components
    ///
    /// Should replace return with in-house iterator struct.
    #[inline]
    pub fn degs_mut(&mut self) -> std::slice::IterMut<D> {
        self.0.iter_mut()
    }

    /// returns the 'total degree' of a multidegree, i.e., the sum of the
    /// individual degree components
    #[inline]
    pub fn total_deg(&self) -> D {
        self.degs().sum()
    }

    /// returns the maximum index for which `self` contains an entry.
    ///
    /// in other words, this is the minimum value `n` for which we would
    /// consider `self` to be an element of the polynomial ring in `n` variables
    #[inline]
    pub fn len(&self) -> I {
        self.0.len()
    }

    pub fn is_succ(&self, other: &MultiDegree) -> bool {
        if self.len() > other.len() {
            return false;
        }

        for (a, b) in self.degs().zip(other.degs()) {
            if a > b {
                return false;
            }
        }

        true
    }

    /// Checked subtraction.
    /// Computes `self - rhs`, returning none if `self[i] < rhs[i]` for any `i`.
    pub fn checked_sub(&self, rhs: &MultiDegree) -> Option<MultiDegree> {
        if rhs.is_zero() {
            return Some(self.clone());
        }
        // `other` is guaranteed nonzero

        if self.is_zero() {
            // `self - other` < 0 so subtraction would fail
            return None;
        }
        // `self` and `other` are guaranteed nonzero

        let mut degs = Vec::with_capacity(self.len());
        let mut zero_cache = 0;

        for pair in self.degs().zip_longest(rhs.degs()) {
            let c = match pair {
                EitherOrBoth::Both(a, b) => match a.cmp(b) {
                    // `a - b < 0`; subtraction has failed
                    Ordering::Less => return None,
                    // `a - b == 0`; increment zero cache
                    Ordering::Equal => {
                        zero_cache += 1;
                        continue;
                    }
                    // `a - b > 0`; carry on with subtraction
                    Ordering::Greater => a - b,
                },
                // only `self` degs remaining, last of which is trusted to be nonzero, so just take the rest
                EitherOrBoth::Left(a) => *a,
                // only `rhs` degs remaining, last of which is trusted to be nonzero, so subtraction has failed
                EitherOrBoth::Right(_) => return None,
            };
            // `c` is guaranteed nonzero

            // catch up with zeros
            degs.append(&mut vec![0; zero_cache]);
            zero_cache = 0;

            degs.push(c);
        }

        Some(MultiDegree(degs))
    }
}

impl Default for MultiDegree {
    #[inline]
    fn default() -> Self {
        Self::zero()
    }
}

#[binop(derefs)]
impl Add for &MultiDegree {
    type Output = MultiDegree;

    fn add(self, rhs: &MultiDegree) -> Self::Output {
        let mut degs = Vec::with_capacity(self.len().max(rhs.len()));
        let mut zero_cache = 0;

        for pair in self.degs().zip_longest(rhs.degs()) {
            let c = match pair {
                EitherOrBoth::Both(0, 0) | EitherOrBoth::Left(0) | EitherOrBoth::Right(0) => {
                    zero_cache += 1;
                    continue;
                }
                EitherOrBoth::Both(a, b) => a + b,
                EitherOrBoth::Left(a) => *a,
                EitherOrBoth::Right(b) => *b,
            };
            // `c` is guaranteed nonzero

            // catch up with zeros
            degs.append(&mut vec![0; zero_cache]);
            zero_cache = 0;

            degs.push(c);
        }

        MultiDegree(degs)
    }
}
