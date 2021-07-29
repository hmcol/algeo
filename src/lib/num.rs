use std::ops::{Add, Sub, Mul, Div};


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


// field -----------------------------------------------------------------------

pub trait Field:
    Sized
    + Add<Self, Output = Self>
    + Sub<Self, Output = Self>
    + Zero
    + Mul<Self, Output = Self>
    + Div<Self, Output = Self>
    + One
{

}

macro_rules! field_impl {
    ($($t:ty)*) => ($(
        impl Field for $t {
            
        }
    )*)
}

field_impl! { f32 f64 }

// pow -------------------------------------------------------------------------
