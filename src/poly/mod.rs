/// Struct and trait implementations relating to the polynomial ring $F[x_1, \dots, x_n]$
pub mod ring;
/// Computation related to finding Gröbner bases.
pub mod comp;
/// Monomial orders.
pub mod ord;



#[doc(inline)]
pub use ring::{u, v, w, x, y, z, MDeg, Polynomial, Term};