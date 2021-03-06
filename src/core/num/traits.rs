use std::cmp::Ordering;
use std::iter::{Product, Sum};
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use super::int::Integer;
use super::rat::Rational;

// zero ------------------------------------------------------------------------

pub trait Zero: Sized + Add<Self, Output = Self> {
    /// if this doesn't work use function, like joseph said
    const ZERO: Self;
    fn is_zero(&self) -> bool;
}

macro_rules! impl_zero_for_primitives {
    ($($t:ty)*) => {
        $(
            impl Zero for $t {
                const ZERO: Self = 0 as Self;

                #[inline]
                fn is_zero(&self) -> bool {
                    *self == Self::ZERO
                }
            }
        )*
    }
}

impl_zero_for_primitives! { usize u8 u16 u32 u64 u128 isize i8 i16 i32 i64 i128 f32 f64 }

impl Zero for Integer {
    const ZERO: Self = Integer(0);

    #[inline]
    fn is_zero(&self) -> bool {
        self.0 == 0
    }
}

impl Zero for Rational {
    const ZERO: Self = Rational {
        numer: Integer::ZERO,
        denom: Integer::ONE,
    };

    #[inline]
    fn is_zero(&self) -> bool {
        self.numer.is_zero()
    }
}

// one -------------------------------------------------------------------------

pub trait One: Sized + Mul {
    /// if this doesn't work use function, like joseph said
    const ONE: Self;
    fn is_one(&self) -> bool;
}

macro_rules! impl_one_for_primitives {
    ($($t:ty)*) => {
        $(
            impl One for $t {
                const ONE: Self = 1 as Self;

                #[inline]
                fn is_one(&self) -> bool {
                    *self == Self::ONE
                }
            }
        )*
    }
}

impl_one_for_primitives! { usize u8 u16 u32 u64 u128 isize i8 i16 i32 i64 i128 f32 f64 }

impl One for Integer {
    const ONE: Self = Integer(1);

    #[inline]
    fn is_one(&self) -> bool {
        self.0 == 1
    }
}

impl One for Rational {
    const ONE: Self = Rational {
        numer: Integer::ONE,
        denom: Integer::ONE,
    };

    #[inline]
    fn is_one(&self) -> bool {
        // due to how `Rational` is implemented, being one always means that numer and denom are both one, so the following could be used and is sort of more semantic
        // self.numer.is_one() && self.denom.is_one()

        // however, doing only one comparison is quicker and more (unnecessarily) general
        // i.e., this would work even if rationals were not simplest form sanitized
        self.numer == self.denom
    }
}
// field -----------------------------------------------------------------------

macro_rules! auto_trait {
    ($Trait:ident: $sup1:path $(, $sup:path)*) => {
        pub trait $Trait: $sup1 $(+ $sup)* {}
        impl<T: $sup1 $(+ $sup)*> $Trait for T {}
    };
}

auto_trait! { OverloadAddition: Add<Self, Output = Self>, AddAssign<Self>, Sum }
auto_trait! { OverloadSubtraction: Sized, Sub<Self, Output = Self>, SubAssign<Self>, Neg<Output = Self> }
auto_trait! { OverloadMultiplication: Mul<Self, Output = Self>, MulAssign<Self>, Product }
auto_trait! { OverloadDivision: Sized, Div<Self, Output = Self>, DivAssign<Self> }

pub trait Field:
    Sized
    + Copy
    + std::fmt::Debug
    + std::fmt::Display
    + PartialEq
    + OverloadAddition
    + OverloadSubtraction
    + Zero
    + OverloadMultiplication
    + OverloadDivision
    + One
{
    fn powi32(&self, p: i32) -> Self;
    /// multiplicative inverse
    fn inv(self) -> Self {
        Self::ONE / self
    }
}

macro_rules! field_impl {
    ($($t:ty)*) => {
        $(
            impl Field for $t {
                fn powi32(&self, exp: i32) -> Self {
                    self.powi(exp)
                }
                fn inv(self) -> Self {
                    self.recip()
                }
            }
        )*
    }
}

field_impl! { f32 f64 Rational }

// Epsilon Equal ---------------------------------------------------------------

pub trait EpsilonEquality {
    fn epsilon_equals(&self, other: &Self) -> bool;
}

macro_rules! epsilon_equality_impl {
    ($($t:ty)*) => {
        $(
            impl EpsilonEquality for $t {
                fn epsilon_equals(&self, other: &Self) -> bool {
                    <$t>::abs(self - other) < <$t>::EPSILON
                }
            }
        )*
    }
}

epsilon_equality_impl! { f32 f64 }

impl EpsilonEquality for Rational {
    fn epsilon_equals(&self, other: &Self) -> bool {
        self == other
    }
}

// Numerical Stability Norm ----------------------------------------------------

/// Partial order to indicate which element has greater stability.
/// More precisely, it is used to determine what number is best to divide by.
/// For lu decomposition (or just Gaussian elimination in general), it is
/// always best to divide by large floats, since this will lead to the least
/// rounding/ float point problems.
///
/// Needs to satisfy 0<= every number, (and every number !<= 0) in order
/// to avoid divide by 0.
pub trait StabilityCmp {
    fn stability_cmp(&self, other: &Self) -> Option<Ordering>;
}

macro_rules! stability_cmp_impl {
    ($($t:ty)*) => {
        $(
            impl StabilityCmp for $t {
                fn stability_cmp(&self, other: &Self) -> Option<Ordering> {
                    self.abs().partial_cmp(&other.abs())
                }
            }
        )*
    }
}

stability_cmp_impl! { f32 f64 }

/// stability cmp for Rational. First priority is making it so that 0 is leq than
/// everything in order for 0 to never be chosen as max value to divide by
/// (that would do division by 0 and get NaN). Otherwise, the ordering chooses
/// a value based on largest sum of abs of numer and denom, in order to have
/// avoid very large denominators/numerators.
impl StabilityCmp for Rational {
    fn stability_cmp(&self, other: &Self) -> Option<Ordering> {
        if *self == Rational::ZERO {
            Some(Ordering::Less)
        } else if *other == Rational::ZERO {
            Some(Ordering::Greater)
        } else {
            (other.numer.abs() + other.denom.abs())
                .partial_cmp(&(self.numer.abs() + self.denom.abs()))
        }
    }
}

impl PartialOrd for Rational {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Rational {
    fn cmp(&self, other: &Self) -> Ordering {
        (self.numer * other.denom).cmp(&(other.numer * self.denom))
    }
}
