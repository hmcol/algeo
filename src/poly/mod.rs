#[macro_use]
pub mod macros;

/// Computation related to finding Gröbner bases.
pub mod comp;
/// Monomial orders.
pub mod ord;
/// Struct and trait implementations relating to the polynomial ring $F[x_1, \dots, x_n]$
pub mod elts;

pub mod mdeg;
mod display;


/* #[doc(inline)]
pub use ring::{u, v, w, x, y, z, Polynomial, Term}; */
