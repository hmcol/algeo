use std::ops::{Add, Sub, Mul, Div};
mod lib;
use lib::poly::*;


fn main() {
    let a = MultiDegree::from_slice(&[2, 3, 4]);
    let b = MultiDegree::from_slice(&[0, 1, 2, 0, 3]);


    println!("{:?}", &a + &b);

    println!("{:?}", a);
}
