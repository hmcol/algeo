use super::super::num::{Field, StabilityCmp, EpsilonEquality};
use super::mat::Mat;
use super::util::get_max_index;


pub struct LUDecomposition<F: Field> {
	pub p : Mat<F>,
	pub l : Mat<F>,
	pub u : Mat<F>
}

pub struct RREFDecomposition<F: Field> {
	pub p : Mat<F>,
	pub rref : Mat<F>
}


impl<F: Field+StabilityCmp+EpsilonEquality> Mat<F> {

	/// Returns LU Decomposition of given matrix. This method does not make
	/// any assumptions on the dimension of the matrix (i.e. we allow
	/// rectangular matrices), though it does require the field to implement 
	/// StabilityCmp in order to do Gaussian elimination under the hood.
	/// 
	/// The LU decomposition satisfies more requirements than that typically
	/// seen in Fortran implementations. The `u` matrix is actually in row
	/// echelon form, which will be convenient for calculating rank as well as
	/// rref (which will be used to compute a basis for the kernel)
	pub fn lu(&self) -> LUDecomposition<F> {
		// Matrix that is eventually made into row echelon form. It
		// will become upper triangular.
		let mut u = self.clone();
		let n = u.rows();

		// elem_prod is product of inverses of elementary matrices applied
		// to `mat`, perm is permutation, such that
		// `perm * self=elem_prod * mat`
		// or equivalently
		// `self = perm^T * elem_prod * mat`
		let mut p_vec : Vec<usize> = (0..n).collect();
		let mut l : Mat<F> = Mat::identity(n);

		// row that hasn't been cleared yet
		let mut r = 0;

		// iterate through all columns
		for c in 0..u.cols() {
			// if cleared out all the rows before we've gone through the
			// columns, just finish here
			if r >= u.rows() {
				break;
			}

			if let Some((stable_index, scalar)) = u.get_stable_column_index_under_row(c, r) {
				// permute
				p_vec.swap(r, stable_index);
				u.permute(r, stable_index);
				l.permute_under_diagonal(r, stable_index);

				// scale l and u correspondingly
				l[(r,r)] = scalar;
				u.scale(r, F::ONE/scalar);

				// then, we do replacement
				for r2 in (r+1)..u.rows() {
					let scalar2 = u[(r2,c)];
					u.replace(r, r2, F::ZERO-scalar2);
					l[(r2, r)] = scalar2;
				}
				
				// go on to the next row, since after doing replacement we should
				// only look beyond in order to get row echelon form.
				r += 1;
			}
		}
		LUDecomposition{
			p: Mat::permutation_from_vec(p_vec),
			l: l,
			u: u
		}
	}

	fn get_stable_column_index_under_row(&self, c: usize, r: usize)-> Option<(usize, F)>
	{
		let stable_index = get_max_index(
			self.get_col_iter(c).skip(r),
			|x,y| x.stability_cmp(y)
		)+r;

		let scalar = self[(stable_index, c)];
		if scalar.epsilon_equals(&F::ZERO) {
			None
		} else {
			Some((stable_index, scalar))
		}
	}

	// assumes src < target
	fn permute_under_diagonal(&mut self, src: usize, target: usize) {
		for c in 0..src {
			let temp = self[(src, c)];
			self[(src, c)] = self[(target, c)];
			self[(target, c)] = temp;
		}
	}

	// pub fn reduced_row_echelon(&self) -> RREFDecomposition {
	// 	let lu = self.lu();

	// 	let mut mat = lu.l;
	// 	let n = mat.rows();

	// 	let mut p = lu.p;

	// 	for r in (1..n).rev() {
	// 		let c_option = temp.index_of_first_nonzero_entry(r);
	// 		if let Some(c) = c_option {
	// 			for r2 in 0..r {
	// 				temp = &Mat::replacement(n, r, r2, F::ZERO-*temp.get_unchecked(r2,c)) * &temp;
	// 			}
	// 		}
	// 	}

	// 	RREFDecomposition {
	// 		p: p,
	// 		rref: mat
	// 	}
	// }
}

#[cfg(test)]
mod tests {
	use super::Mat;
	use super::super::super::num::EpsilonEquality;
	use super::super::util::num_to_seq;
	
	// #[test]
	// fn test_lu_1(){
	// 	let mat : Mat<f32> = Mat::new(5,5,
	// 		vec![
	// 			5.0, 5.0, 5.0,  5.0,  1.0,
	// 			4.0, 5.0, 6.0,  7.0,  1.0,
	// 			8.0, 9.0, 10.0, 11.0, 1.0,
	// 			0.0, 9.0, 2.0,  3.0,  1.0,
	// 			7.0, 6.0, 9.0,  4.0,  1.0
	// 		]
	// 	);
	// 	let lu = mat.lu();

	// 	assert_eq!(mat, &lu.p.transpose() * &(&lu.l * &lu.u));
	// }

	#[test]
	fn bad_case() {
		let mut mat : Mat<f64> = Mat::new(3, 4,
			vec![
				1.0, 2.0, 2.0, 0.0,
				2.0, 0.0, 1.0, 0.0,
				2.0, 1.0, 0.0, 0.0
			]
		);
		let lu = mat.lu();
		println!("{}", lu.u);
		println!("{}", lu.l);
		println!("{}", lu.p);
		assert!(mat.epsilon_equals(&(&lu.p.transpose() * &(&lu.l * &lu.u))));
	}

	#[test]
	fn test_lu() {

		let n: usize = 3;
		let m: usize = 4;
		let values = vec![0.0, 1.0, 2.0];

		// iterate through all matrices with only entries in `values`
		for v in 0..(values.len().pow((n*m) as u32)) {
			let mut entries: Vec<f32> = num_to_seq(v, values.len())
				.map(|i| values[i]).collect();
			
			for _i in (entries.len())..(n*m) {
				entries.push(0.0);
			}

			let mat : Mat<f32> = Mat::new(n, m, entries);
			let lu = mat.lu();
			let prod = &lu.p.transpose()*&(&lu.l*&lu.u);

			assert!(lu.l.is_lower_triangular(),
				"L={} is not lower triangular, original A={}",
				lu.l, mat
			);
			assert!(lu.u.is_upper_triangular(),
				"U={} is not upper triangular, original A={}",
				lu.u, mat
			);
			assert!(mat.epsilon_equals(&prod),
				"failed on {} == {}, epsilon error is {}",
				mat, prod, (&mat + &(&prod*(-1.0))).frobenius_norm()
			);
		}
	}
}