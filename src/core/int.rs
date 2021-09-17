use std::{collections::BTreeMap, ops::{Add, Div, Mul, Neg, Rem, Sub}};

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Default, Debug)]
struct Integer(i64);

macro_rules! pull_unop {
    ($($func:ident)*) => {
        $(
            #[inline]
            fn $func(self) -> Integer {
                Integer(self.0.$func())
            }
        )*
    };
}

macro_rules! pull_check {
    ($($func:ident)*) => {
        $(
            #[inline]
            fn $func(self) -> bool {
                self.0.$func()
            }
        )*
    };
}

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
    pull_unop! { abs signum }
    pull_check! { is_positive is_negative }

    #[inline]
    fn is_zero(self) -> bool {
        self.0 == 0
    }

    #[inline]
    fn is_one(self) -> bool {
        self.0 == 1
    }


    fn gcd(a: Integer, b: Integer) -> Integer {
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

    pull_unop! { neg }
}


#[cfg(test)]
mod tests {
    use super::Integer;

    fn integer() {}
}
