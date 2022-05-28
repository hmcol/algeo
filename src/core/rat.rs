#![allow(unused)]


use std::{
    iter::{Product, Sum},
    ops::{Add, Div, Mul, Neg, Sub},
};

use super::int::{gcd, lcm, Integer};

#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub struct Rational {
    pub numer: Integer,
    pub denom: Integer,
}

impl Rational {
    const ZERO: Rational = Rational::new_unchecked(0, 1);
    const ONE: Rational = Rational::new_unchecked(1, 1);

    /// returns rational `numer`/`denom` without simplifying
    ///
    /// not for external use; rationals must always be created in simplest form
    #[inline]
    const fn new_unchecked(numer: Integer, denom: Integer) -> Self {
        Rational { numer, denom }
    }

    pub fn new(numer: Integer, denom: Integer) -> Self {
        if denom == 0 {
            panic!("denominator cannot be zero")
        };

        let gcd = gcd(numer, denom);

        Self::new_unchecked(denom.signum() * numer / gcd, denom.abs() / gcd)
    }

    /// returns rational `int/1`
    #[inline]
    pub fn new_int(int: Integer) -> Self {
        Self::new_unchecked(int, 1)
    }

    // pub fn pow(&self, exp: Integer) -> Self {
    //     if exp.is_neg() {
    //         Self::new_unchecked(self.denom.pow(-exp), self.numer.pow(-exp))
    //     } else {
    //         Self::new_unchecked(self.numer.pow(exp), self.denom.pow(exp))
    //     }
    // }

    // /// for compatibility with Field trait
    // pub fn powi(&self, exp: i32) -> Self {
    //     self.pow(exp.into())
    // }

    pub fn signum(&self) -> Integer {
        self.numer.signum()
    }

    pub fn recip(&self) -> Self {
        Self::new_unchecked(self.signum() * self.denom, self.signum() * self.numer)
    }
}

impl From<Integer> for Rational {
    fn from(int: Integer) -> Self {
        Self::new_int(int)
    }
}

impl Add for Rational {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        let numer = self.numer * other.denom + other.numer * self.denom;

        if numer == 0 {
            return Rational::ZERO;
        }

        Rational::new_unchecked(numer, self.denom * other.denom)
    }
}

impl Sum for Rational {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.reduce(Add::add).unwrap_or(Rational::ZERO)
    }
}

impl Sub for Rational {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        let numer = self.numer * other.denom - other.numer * self.denom;

        if numer == 0 {
            return Rational::ZERO;
        }

        Rational::new_unchecked(numer, self.denom * other.denom)
    }
}

impl Neg for Rational {
    type Output = Self;

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
        iter.reduce(Mul::mul).unwrap_or(Rational::ONE)
    }
}

/// uhhh, who should handle dividing by zero?
impl Div for Rational {
    type Output = Self;

    fn div(self, other: Self) -> Self {
        Rational::new(self.numer * other.denom, self.denom * other.numer)
    }
}
