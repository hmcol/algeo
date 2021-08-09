use super::super::num::{Field, StabilityCmp, EpsilonEquality};
use super::mat::Mat;
use super::util::get_max_index;

pub struct LUDecomposition<F: Field> {
	pub p : Mat<F>,
	pub l : Mat<F>,
	pub u : Mat<F>
}

impl<F: Field+StabilityCmp+EpsilonEquality> Mat<F> {

	/// NOT TESTED AT ALL, PROBABLY LOTS OF BUGS.
	/// 
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
		let mut mat = self.clone();
		let n = mat.rows();

		// elem_prod is product of inverses of elementary matrices applied
		// to `mat`, perm is permutation, such that
		// `perm * self=elem_prod * mat`
		// or equivalently
		// `self = perm^T * elem_prod * mat`
		let mut perm : Vec<usize> = (0..n).collect();
		let mut elem_prod : Mat<F> = Mat::identity(n);

		// row that hasn't been cleared yet
		let mut r = 0;

		// iterate through all columns
		for c in 0..mat.cols() {
			// if cleared out all the rows before we've gone through the
			// columns, just finish here
			if r >= mat.rows() {
				break;
			}

			// get index of most stable (usually largest) entry in column after
			// the current row.
			let stable_index = get_max_index(
				mat.get_col_iter(c).skip(r),
				|x,y| x.stability_cmp(y)
			)+r;

			// if "most stable" index is zero, then no replacement is needed
			// in this column, so we skip to the next.
			//
			// TODO: switch from `scalar == F::ZERO` to something like
			// `scalar-F::ZERO < eps`. Somehow this should be integrated in
			// StabilityCmp.
			let scalar = mat[(stable_index, c)];
			if scalar.epsilon_equals(&F::ZERO) {
				continue;
			}

			// otherwise, we permute both our matrix and our permutation
			// to get the stable_index to be in the appropriate row.
			perm.swap(r, stable_index);
			mat.permute(r, stable_index);

			// permute L matrix, but only entries under the diagonal
			if r>0 {
				for c2 in 0..r {
					let temp = elem_prod[(r, c2)];
					elem_prod[(r, c2)] = elem_prod[(stable_index, c2)];
					elem_prod[(stable_index, c2)] = temp;
				}
			}

			// scale l and u correspondingly
			elem_prod[(r,r)] = scalar;
			mat.scale(r, F::ONE/scalar);

			// then, we do replacement
			for r2 in (r+1)..mat.rows() {
				let scalar2 = mat[(r2,c)];
				mat.replace(r, r2, F::ZERO-scalar2);
				elem_prod[(r2, r)] = scalar2;
			}
			
			// go on to the next row, since after doing replacement we should
			// only look beyond in order to get row echelon form.
			r += 1;
		}
		LUDecomposition{
			p: Mat::permutation_from_vec(perm),
			l: elem_prod,
			u: mat
		}
	}
}

#[cfg(test)]
mod tests {
	use super::Mat;
	use super::super::super::num::EpsilonEquality;
	use super::super::util::num_to_seq;
	
	#[test]
	fn test_lu_1(){
		let mat : Mat<f32> = Mat::new(5,5,
			vec![
				5.0, 5.0, 5.0,  5.0,  1.0,
				4.0, 5.0, 6.0,  7.0,  1.0,
				8.0, 9.0, 10.0, 11.0, 1.0,
				0.0, 9.0, 2.0,  3.0,  1.0,
				7.0, 6.0, 9.0,  4.0,  1.0
			]
		);
		let lu = mat.lu();

		assert_eq!(mat, &lu.p.transpose() * &(&lu.l * &lu.u));
	}

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
				"L={} is not lower triangular",
				lu.l
			);
			assert!(lu.u.is_upper_triangular(),
				"U={} is not upper triangular",
				lu.u
			);
			assert!(mat.epsilon_equals(&prod),
				"failed on {} == {}",
				mat, prod
			);
		}
	}
}