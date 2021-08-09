
mod lib;
use lib::poly::*;
use lib::linalg::mat::Mat;
use lib::linalg::row_echelon::*;


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



    let mat : Mat<f32> = Mat::new(3, 4,
        vec![
            0.0, 0.0, 0.0, 1.0,
            0.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 0.0
        ]
    );
    let lu = mat.lu();
    println!("{}", lu.u);
    println!("{}", lu.l);
    println!("{}", lu.p);

	//println!("{:?}", mat.reduced_row_echelon().0);
	//println!("{:?}", mat.compute_kernel());

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
            let grlex = cmp_grlex(a, b);
            let grevlex = cmp_grevlex(a, b);

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