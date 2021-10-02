use std::{
    iter::{Product, Sum},
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

use super::{
    Integer,
    traits::{One, Zero},
};

#[derive(Copy, Clone, PartialEq, Eq, Debug, Hash)]
pub struct Rational {
    pub numer: Integer,
    pub denom: Integer,
}

impl Rational {
    /// returns rational `numer`/`denom` without simplifying
    ///
    /// not for external use; rationals must always be created in simplest form
    #[inline]
    const fn new_unchecked(numer: Integer, denom: Integer) -> Self {
        Rational { numer, denom }
    }

    pub fn new(numer: Integer, denom: Integer) -> Self {
        if denom.is_zero() {
            panic!("denominator cannot be zero")
        };

        let gcd = Integer::gcd(numer, denom);

        Self::new_unchecked(denom.sgn() * numer / gcd, denom.abs() / gcd)
    }

    /// returns rational `numer`/`denom` in simplest form
    #[inline]
    pub fn new_i64(numer: i64, denom: i64) -> Self {
        Self::new(Integer::from(numer), Integer::from(denom))
    }

    /// returns rational `int/1`
    #[inline]
    pub fn new_int(int: Integer) -> Self {
        Self::new_unchecked(int, Integer::ONE)
    }

    pub fn pow(&self, exp: Integer) -> Self {
        if exp.is_neg() {
            Self::new_unchecked(self.denom.pow(-exp), self.numer.pow(-exp))
        } else {
            Self::new_unchecked(self.numer.pow(exp), self.denom.pow(exp))
        }
    }

    /// for compatibility with Field trait
    pub fn powi(&self, exp: i32) -> Self {
        self.pow(Integer(exp.into()))
    }

    #[inline]
    pub fn sgn(&self) -> Integer {
        self.numer.sgn()
    }

    pub fn recip(&self) -> Self {
        Self::new_unchecked(self.sgn() * self.denom, self.sgn() * self.numer)
    }
}

impl From<Integer> for Rational {
    #[inline]
    fn from(int: Integer) -> Self {
        Self::new_int(int)
    }
}

impl Add for Rational {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        let numer = self.numer * other.denom + other.numer * self.denom;

        if numer.is_zero() {
            Rational::ZERO
        } else {
            Rational::new_unchecked(numer, self.denom * other.denom)
        }
    }
}

impl Sum for Rational {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.reduce(Add::add).unwrap_or_else(|| Rational::ZERO)
    }
}

impl Sub for Rational {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        let numer = self.numer * other.denom - other.numer * self.denom;

        if numer.is_zero() {
            Rational::ZERO
        } else {
            Rational::new_unchecked(numer, self.denom * other.denom)
        }
    }
}

impl Neg for Rational {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self::Output {
        Rational::new_unchecked(-self.numer, self.denom)
    }
}

impl Mul for Rational {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        Rational::new(self.numer * other.numer, self.denom * other.denom)
    }
}

impl Product for Rational {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.reduce(Mul::mul).unwrap_or_else(|| Rational::ONE)
    }
}

impl Div for Rational {
    type Output = Self;

    fn div(self, other: Self) -> Self {
        Rational::new(self.numer * other.denom, self.denom * other.numer)
    }
}

/// implements the given operation-assignment traits for rationals in the most
/// obvious way
///
/// for example if the operation in question is star, `⋆`, then `a ⋆= b`
/// is the same as `a = a ⋆ b`.
macro_rules! impl_op_assign_for_rational_simple {
    ($trait_name:ty, $op_ass_func:ident, $op_func:ident) => {
        impl $trait_name for Rational {
            #[inline]
            fn $op_ass_func(&mut self, other: Self) {
                *self = self.$op_func(other);
            }
        }
    };
}

impl_op_assign_for_rational_simple! { AddAssign, add_assign, add }
impl_op_assign_for_rational_simple! { SubAssign, sub_assign, sub }
impl_op_assign_for_rational_simple! { MulAssign, mul_assign, mul }
impl_op_assign_for_rational_simple! { DivAssign, div_assign, div }

impl std::fmt::Display for Rational {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        if self.numer.is_zero() {
            write!(f, "0")
        } else if self.denom.is_one() {
            write!(f, "{}", self.numer)
        } else {
            write!(f, "{}/{}", self.numer, self.denom)

            // bonus latex mode
            // write!(f, "\\frac{{{}}}{{{}}}", self.numer, self.denom)
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::core::num::int::Integer;
    use crate::core::num::traits::One;
    use crate::core::num::traits::Zero;

    use super::Rational;

    #[test]
    fn test_rational() {
        let int = |i: i64| Integer::from(i);

        assert_eq!(
            Rational::new_i64(0, 1),
            Rational {
                numer: Integer::ZERO,
                denom: Integer::ONE
            }
        );
        assert_eq!(
            Rational::new_i64(1, 1),
            Rational {
                numer: Integer::ONE,
                denom: Integer::ONE
            }
        );
        assert_eq!(
            Rational::new_i64(2, 2),
            Rational {
                numer: Integer::ONE,
                denom: Integer::ONE
            }
        );
        assert_eq!(
            Rational::new_i64(1, 2) + Rational::new_i64(1, 3),
            Rational {
                numer: int(5),
                denom: int(6)
            }
        );
        assert_eq!(
            Rational::new_i64(1, 2) - Rational::new_i64(1, 3),
            Rational {
                numer: int(1),
                denom: int(6)
            }
        );
        assert_eq!(
            Rational::new_i64(1, 2) * Rational::new_i64(1, 3),
            Rational {
                numer: int(1),
                denom: int(6)
            }
        );
        assert_eq!(
            Rational::new_i64(1, 2) / Rational::new_i64(1, 3),
            Rational {
                numer: int(3),
                denom: int(2)
            }
        );
    }
}
