use std::ops::{Add, AddAssign, Sub, SubAssign, Mul, MulAssign, Div, DivAssign};

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
    + Add<Self, Output = Self> + AddAssign<Self>
    + Sub<Self, Output = Self> + SubAssign<Self>
    + Zero
    + Mul<Self, Output = Self> + MulAssign<Self>
    + Div<Self, Output = Self> + DivAssign<Self>
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

