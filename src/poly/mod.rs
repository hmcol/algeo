/// elements of the polynomial ring `F[X]`
mod ring;
mod ord;

#[doc(inline)]
pub use ring::{u, v, w, x, y, z, Const, MDeg, Polynomial, Term};
pub use ord::*;