use algeo::{
    core::num::Rational,
    poly::{comp::Computer, elts::Polynomial, ord},
};

type Poly = Polynomial<Rational>;
type Comp = Computer<Rational, ord::GRevLex>;

// use algeo::linalg::mat::Mat;
//use algeo::linalg::row_echelon::*;

macro_rules! fn_vars {
    (@idx x) => { 0 };
    (@idx y) => { 1 };
    (@idx z) => { 2 };
    (@idx w) => { 3 };
    (@idx u) => { 4 };
    (@idx v) => { 5 };
    ($($var:ident)*) => {
        $(
            fn $var(deg: u8) -> Poly {
                Poly::var(fn_vars!(@idx $var), deg)
            }
        )*
    };
}

fn main() {
    fn c(coef: i64) -> Poly {
        Poly::from(Rational::new_i64(coef, 1))
    }
    fn_vars! { x y z }

    let mut gens: Vec<Poly> = vec![
        x(4) - y(4) + z(3) - c(1),
        x(3) + y(2) + z(2) - c(1),
    ];

    println!("INITIAL");
    for (i, g) in gens.iter().enumerate() {
        println!("g[{}] = {}", i, g);
    }

    Comp::buchberger_algorithm(&mut gens);

    println!("\nFINAL");
    for (i, g) in gens.iter().enumerate() {
        println!("g[{}] = {}", i, g);
    }

}
