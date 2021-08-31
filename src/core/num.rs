use std::ops::{Add, Sub, Mul, Div};
use std::cmp::Ordering;

use super::frac::Frac;


// zero ------------------------------------------------------------------------

pub trait Zero: Sized + Add<Self, Output = Self> {
    /// if this doesn't work use function, like joseph said
    const ZERO: Self;
}

macro_rules! int_zero_impl {
    ($($t:ty)*) => ($(
        impl Zero for $t {
            const ZERO: Self = 0 as Self;
        }
    )*)
}

int_zero_impl! { usize u8 u16 u32 u64 u128 isize i8 i16 i32 i64 i128 f32 f64 }

impl Zero for Frac {
    const ZERO: Frac = Frac::zero();
}

// one -------------------------------------------------------------------------

pub trait One: Sized + Mul {
    /// if this doesn't work use function, like joseph said
    const ONE: Self;
}

macro_rules! one_impl {
    ($($t:ty)*) => ($(
        impl One for $t {
            const ONE: Self = 1 as Self;
        }
    )*)
}

one_impl! { usize u8 u16 u32 u64 u128 isize i8 i16 i32 i64 i128 f32 f64 }

impl One for Frac {
    const ONE: Frac = Frac::one();
}

// field -----------------------------------------------------------------------

pub trait Field:
    Sized
    + Copy
    + std::fmt::Debug
    + std::fmt::Display
    + PartialEq
    + Add<Self, Output = Self>
    + Sub<Self, Output = Self>
    + Zero
    + Mul<Self, Output = Self>
    + Div<Self, Output = Self>
    + One
{
    fn powi32(&self, p: i32) -> Self;
}

macro_rules! field_impl {
    ($($t:ty)*) => ($(
        impl Field for $t {
            fn powi32(&self, p: i32) -> Self {
                self.powi(p)
            }
        }
    )*)
}

field_impl! { f32 f64 Frac }

// Epsilon Equal ---------------------------------------------------------------

pub trait EpsilonEquality {
    fn epsilon_equals(&self, other: &Self) -> bool;
}

macro_rules! epsilon_equality_impl {
    ($($t:ty)*) => ($(
        impl EpsilonEquality for $t {
            fn epsilon_equals(&self, other: &Self) -> bool {
                (self-other).abs() < <$t>::EPSILON
            }
        }
    )*)
}

epsilon_equality_impl! { f32 f64 }

impl EpsilonEquality for Frac {
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
    ($($t:ty)*) => ($(
        impl StabilityCmp for $t {
            fn stability_cmp(&self, other: &Self) -> Option<Ordering> {
                self.abs().partial_cmp(&other.abs())
            }
        }
    )*)
}

stability_cmp_impl! { f32 f64 }

/// stability cmp for Frac. First priority is making it so that 0 is leq than
/// everything in order for 0 to never be chosen as max value to divide by
/// (that would do division by 0 and get NaN). Otherwise, the ordering chooses
/// a value based on largest sum of abs of numer and denom, in order to have
/// avoid very large denominators/numerators.
impl StabilityCmp for Frac {
    fn stability_cmp(&self, other: &Self) -> Option<Ordering> {
        if *self == Frac::ZERO {
            Some(Ordering::Less)
        } else if *other == Frac::ZERO {
            Some(Ordering::Greater)
        } else {
            (other.numer.abs()+other.denom.abs())
                .partial_cmp(&(self.numer.abs()+self.denom.abs()))
        }
    }
}

impl PartialOrd for Frac {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Frac {
    fn cmp(&self, other: &Self) -> Ordering {
        (self.numer*other.denom).cmp(&(other.numer*self.denom))
    }
}