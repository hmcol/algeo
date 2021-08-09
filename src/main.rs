
mod lib;
use lib::poly::*;


fn main() {
    let x = Term::monic(MDeg([1, 0]));
    let y = Term::monic(MDeg([0, 1]));
    let neg_y = Term::from_coef_mdeg(-1.0, MDeg([0, 1]));
    
    let mut p = Polynomial::zero();
    let mut q = Polynomial::zero();

    p += x;
    p += y;

    q += x;
    q += neg_y;

    print!("pq = ({})({}) = ", p, q);

    println!("{}", p * q);


}

fn grlex_vs_grevlex() {
    let mut degs = Vec::new();

    const D: i32 = 5;

    for x in 0..D {
        for y in 0..(D - x) {
            for z in 0..(D - x - y) {
                degs.push(MDeg([x, y, z]));
            }
        }
    }
    
    for &a in &degs {
        for &b in &degs {
            let grlex = mon_ord::grlex(a, b);
            let grevlex = mon_ord::grevlex(a, b);

            if grlex != grevlex {
                println!("cmp {} {} => grlex: {:<10} grevlex: {:<10}", 
                    a, b, 
                    format!("{:?}", grlex),
                    format!("{:?}", grevlex)
                );
            }
        }
    }
}