use std::ops::{Add, Sub, Mul, Div};


/// integer type for fractions
pub type I = i64;

#[derive(Copy, Clone, PartialEq, Eq, Debug, Hash)]
pub struct Frac {
    numer: I,
    denom: I,
}

impl Frac {
    /// returns fraction numer/denom without simplifying
    /// 
    /// not for external use; fractions must always be created in simplest form
    const fn new_unchecked(numer: I, denom: I) -> Self {
        Frac {
            numer,
            denom,
        }
    }

    /// returns the zero element: 0/1
    pub const fn zero() -> Self {
        Frac::new_unchecked(0, 1)
    }

    /// returns the one element: 1/1
    pub const fn one() -> Self {
        Frac::new_unchecked(1, 1)
    }

    /// returns fraction numer/denom in simplest form
    pub fn new(numer: I, denom: I) -> Self {
        if denom == 0 { panic!("denominator cannot be zero") };

        let gcd = gcd(numer, denom);

        Frac::new_unchecked(
            numer * denom.signum() / gcd ,
            denom.abs() / gcd,
        )
    }

    /// for compatibility with Field trait
    pub fn powi(&self, n: i32) -> Self {
        Frac::new(
            self.numer.pow(n as u32),
            self.denom.pow(n as u32)
        )
    }
}

impl Add for Frac {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        Frac::new(
            self.numer * rhs.denom + rhs.numer * self.denom,
            self.denom * rhs.denom,
        )
    }
}

impl Sub for Frac {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        Frac::new(
            self.numer * rhs.denom - rhs.numer * self.denom,
            self.denom * rhs.denom,
        )
    }
}

impl Mul for Frac {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        Frac::new(
            self.numer * rhs.numer,
            self.denom * rhs.denom
        )
    }
}

impl Div for Frac {
    type Output = Self;

    fn div(self, rhs: Self) -> Self {
        Frac::new(
            self.numer * rhs.denom,
            self.denom * rhs.numer,
        )
    }
}

impl std::fmt::Display for Frac {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}/{}", self.numer, self.denom)
    }
}

/// computes greatest common divisor using Euclidean algorithm
fn gcd(a: i64, b: i64) -> i64 {
    let mut a = a;
    let mut b = b;
    let mut t: i64;

    while b != 0 {
        t = b;
        b = a % b;
        a = t;
    }

    // return 
    a.abs()
}


#[cfg(test)]
mod tests {
    use super::gcd;
    use super::Frac;

    #[test]
    fn test_gcd() {
        assert_eq!(gcd(-10,-15), 5);
        assert_eq!(gcd(-5,-10), 5);
        assert_eq!(gcd(-0,-1), 1);
        assert_eq!(gcd(-5,-6), 1);
        assert_eq!(gcd(-80,-20), 20);
    }

    #[test]
    fn test_frac() {
        assert_eq!(Frac::new(0,1), Frac{ numer: 0, denom: 1});
        assert_eq!(Frac::new(1,1), Frac{ numer: 1, denom: 1});
        assert_eq!(Frac::new(2,2), Frac{ numer: 1, denom: 1});
        assert_eq!(Frac::new(1,2)+Frac::new(1,3), Frac{ numer: 5, denom: 6});
        assert_eq!(Frac::new(1,2)-Frac::new(1,3), Frac{ numer: 1, denom: 6});
        assert_eq!(Frac::new(1,2)*Frac::new(1,3), Frac{ numer: 1, denom: 6});
        assert_eq!(Frac::new(1,2)/Frac::new(1,3), Frac{ numer: 3, denom: 2});
    }
}