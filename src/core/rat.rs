#![allow(unused)]

use std::{
    iter::{Product, Sum},
    ops::{Add, Div, Mul, Neg, Sub},
};

use super::int::{gcd, lcm, Integer};

#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub struct Rational {
    numer: Integer,
    denom: Integer,
}

impl Rational {
    pub const ZERO: Rational = Rational { numer: 0, denom: 1 };
    pub const ONE: Rational = Rational { numer: 1, denom: 1 };

    pub fn new(numer: Integer, denom: Integer) -> Self {
        if denom == 0 {
            panic!("denominator cannot be zero")
        };

        let gcd = gcd(numer, denom);

        Rational {
            numer: denom.signum() * numer / gcd,
            denom: denom.abs() / gcd,
        }
    }

    /// returns rational `int/1`
    #[inline]
    pub fn new_int(int: Integer) -> Self {
        Rational {
            numer: int,
            denom: 1,
        }
    }

    pub fn is_zero(&self) -> bool {
        self.numer == 0
    }

    pub fn signum(&self) -> Integer {
        self.numer.signum()
    }

    pub fn recip(&self) -> Option<Self> {
        if self.is_zero() {
            return None;
        }

        Some(Rational {
            numer: self.signum() * self.denom,
            denom: self.signum() * self.numer,
        })
    }
}

impl From<Integer> for Rational {
    fn from(int: Integer) -> Self {
        Rational {
            numer: int,
            denom: 1,
        }
    }
}

impl Add for Rational {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        let numer = self.numer * other.denom + other.numer * self.denom;

        if numer == 0 {
            return Rational::ZERO;
        }

        Rational {
            numer,
            denom: self.denom * other.denom,
        }
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

        Rational {
            numer,
            denom: self.denom * other.denom,
        }
    }
}

impl Neg for Rational {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Rational {
            numer: -self.numer,
            denom: self.denom,
        }
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
