use std::ops::{Add, Sub, Mul, Div};

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



// Numerical Stability Norm ----------------------------------------------------

/// Assigns a value in f64 that indicates how "big" a number is. More precisely,
/// it is used to determine what number is best to divide by. For lu
/// decomposition (or just Gaussian elimination in general), it is always best
/// to divide by large floats, since this will lead to the least rounding/
/// float point problems.
pub trait NumStabilityNormed {
    fn num_stability_norm(&self) -> f64;
}

macro_rules! num_stability_norm_impl {
    ($($t:ty)*) => ($(
        impl NumStabilityNormed for $t {
            fn num_stability_norm(&self) -> f64 {
                *self as f64
            }
        }
    )*)
}

num_stability_norm_impl! { f32 f64 }

/// the max stability norm of a vec of fractions will be the one with the
/// one with the minimum sum of numerator+denominator.
impl NumStabilityNormed for Frac {
    fn num_stability_norm(&self) -> f64 {
        (-self.numer-self.denom) as f64
    }
}