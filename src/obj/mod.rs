#![allow(unused)]

mod field;
mod ring;

pub use ring::Ring;
pub use field::Field;

// objects

trait Object {}

pub struct Z;
pub struct ZModNZ {
    modulus: i32,
}
struct ZModPZ {
    modulus: i32,
}

impl ZModPZ {
    /// only returns Z/pZ when p is a prime
    fn new(p: i32) -> Option<ZModPZ> {
        // this should be a check if p is prime
        if true {
            return None;
        }

        Some(ZModPZ { modulus: p })
    }
}

pub struct Q;
