/// elements of the polynomial ring $F[x_1, \dots, x_n]$
pub mod ring_old;
/// general polynomial division
pub mod div;
/// monomial orders
pub mod ord;

mod ring;

#[doc(inline)]
pub use ring::{u, v, w, x, y, z, MDeg, Polynomial, Term};