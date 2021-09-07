//! Polynomial stuff
//! 
//! latex testing
//! 
//! $$ 2^8 $$
//! 
//! $$ \frac{1}{2} $$
//! 
//! $$ \sum_{i = 1}^{n} a_i $$
//! 
//! $$ \mathbb{C}[x_1, \dots, x_n] $$
//! 
//! $$ I = \langle f_1, \dots, f_n \rangle $$


/// Polynomial stuff
/// 
/// latex testing
/// 
/// $$ 2^8 $$
/// 
/// $$ \frac{1}{2} $$
/// 
/// $$ \sum_{i = 1}^{n} a_i $$
/// 
/// $$ \mathbb{C}[x_1, \dots, x_n] $$
/// 
/// $$ I = \langle f_1, \dots, f_n \rangle $$
pub mod ring;
/// general polynomial division
pub mod div;
/// monomial orders
pub mod ord;

#[doc(inline)]
pub use ring::{u, v, w, x, y, z, MDeg, Polynomial, Term};