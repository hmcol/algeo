
use algeo::poly::poly::*;
use algeo::linalg::mat::Mat;
//use algeo::linalg::row_echelon::*;

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
