use std::ops::{Add, Div, Mul, Neg, Rem, Sub};



#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Default, Debug, Hash)]
pub struct Integer(pub i64);


/// possible unary ops ( fn op(self) -> Self ):
/// - abs
/// - signum
/// - factorial
/// - (rounded/floor/ceil) root
/// 
/// possible checks/conditions ( fn is_condition(self) -> bool ):
/// - zero, one
/// - positive, negative
/// - odd, even
/// - irreducible/prime
/// - square, squarefree
/// 
/// possible misc
/// - pow
/// - lcm
impl Integer {
    /// Signum of an integer:
    /// - `-1` if negative
    /// - `0` if zero
    /// - `+1` if positive
    #[inline]
    pub fn sgn(self) -> Self {
        Integer(self.0.signum())
    }

    /// Check whether `self` is positive.
    #[inline]
    pub fn is_pos(self) -> bool {
        self.0.is_positive()
    }

    /// Check whether `self` is negative.
    #[inline]
    pub fn is_neg(self) -> bool {
        self.0.is_negative()
    }

    /// Absolute value of an integer:
    /// - `self` for nonnegative
    /// - `-self` for negative
    #[inline]
    pub fn abs(self) -> Self {
        Integer(self.0.abs())
    }

    /// Computes `self` to the power of `exp`
    /// 
    /// NOTE: `exp` is cast to a `u32` before computation
    #[inline]
    pub fn pow(self, exp: Integer) -> Self {
        Integer(self.0.pow(exp.0 as u32))
    }

    /// Computes the greatest common divisor of `a` and `b`.
    /// 
    /// That is, the largest nonnegative integer which divides both.
    pub fn gcd(a: Self, b: Self) -> Self {
        let mut a = a.0;
        let mut b = b.0;
        let mut t: i64;

        while b != 0 {
            t = b;
            b = a % b;
            a = t;
        }

        a.abs().into()
    }

    /// Computes the least common multiple of `a` and `b`.
    /// 
    /// That is, the smallest nonnegative integer which both divide.
    pub fn lcm(a: Integer, b: Integer) -> Integer {
        a / Self::gcd(a, b) * b
    }
}

impl From<i64> for Integer {
    #[inline]
    fn from(val: i64) -> Self {
        Integer(val)
    }
}



macro_rules! impl_integer_binop {
    ($($Trait:ident, $func:ident);*) => {
        $(
            impl $Trait for Integer {
                type Output = Self;

                #[inline]
                fn $func(self, rhs: Self) -> Self {
                    Integer($Trait::$func(self.0, rhs.0))
                }
            }
        )*
    };
}
impl_integer_binop! { Add, add; Sub, sub; Mul, mul; Div, div; Rem, rem }

impl Neg for Integer {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Integer(self.0.neg())
    }
}


impl std::fmt::Display for Integer {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

#[cfg(test)]
mod tests {
    use super::Integer;

    #[test]
    fn test_gcd() {
        let gcd = |a, b| Integer::gcd(Integer(a), Integer(b)).0;

        assert_eq!(gcd(-10, -15), 5);
        assert_eq!(gcd(-5, -10), 5);
        assert_eq!(gcd(-0, -1), 1);
        assert_eq!(gcd(-5, -6), 1);
        assert_eq!(gcd(-80, -20), 20);
    }

    #[test]
    fn integer() {


    }
}
