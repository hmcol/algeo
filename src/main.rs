
mod lib;
use lib::poly_dyn::*;


fn main() {
    let neg_y = Term::constant(-1.0) * y(1);
    
    let mut p = Polynomial::<f64>::zero();
    let mut q = Polynomial::<f64>::zero();

    p += x(1);
    p += y(1);

    q += x(1);
    q += neg_y;

    println!("pq = ({})({}) = ", p, q);

    println!("{}", p * q);


}