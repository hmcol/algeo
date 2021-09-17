use algeo::poly::elts::Polynomial;

// use algeo::linalg::mat::Mat;
//use algeo::linalg::row_echelon::*;

fn main() {
    let x = |d| Polynomial::var(0, d);
    let y = |d| Polynomial::var(1, d);

    let p: Polynomial<f64> = x(1) + y(1);
    let q: Polynomial<f64> = x(1) - y(1);

    println!("\n({})({}) = {}", p, q, &p * &q);
    // let mat : Mat<f32> = Mat::new(3, 4,
    //     vec![
    //         0.0, 0.0, 0.0, 1.0,
    //         0.0, 0.0, 1.0, 0.0,
    //         0.0, 0.0, 0.0, 0.0
    //     ]
    // );
    // let lu = mat.lu();
    // println!("{}", lu.u);
    // println!("{}", lu.l);
    // println!("{}", lu.p);

    //println!("{:?}", mat.reduced_row_echelon().0);
    //println!("{:?}", mat.compute_kernel());
}
