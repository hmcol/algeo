/// elements of the polynomial ring $F[x_1, \dots, x_n]$
pub mod ring;
/// general polynomial division
pub mod div;
/// monomial orders
pub mod ord;

#[doc(inline)]
pub use ring::{u, v, w, x, y, z, MDeg, Polynomial, Term};