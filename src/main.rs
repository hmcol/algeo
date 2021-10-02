use algeo::{core::num::Rational, poly::elts::Polynomial};

// use algeo::linalg::mat::Mat;
//use algeo::linalg::row_echelon::*;

fn main() {
    let x = |d| Polynomial::<Rational>::var(0, d);
    let y = |d| Polynomial::<Rational>::var(1, d);

    let p = x(1) + y(1);
    let q = x(1) - y(1);

    println!("\n({})({}) = {}", p, q, &p * &q);
}
