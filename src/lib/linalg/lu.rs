use super::super::num::{Field, StabilityCmp};
use super::mat::Mat;
use super::util::get_max_index;

pub struct LUDecomposition<F: Field> {
	pub p : Mat<F>,
	pub l : Mat<F>,
	pub u : Mat<F>
}

impl<F: Field+StabilityCmp> Mat<F> {

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
			let scalar = mat[(stable_index, c)];
			if scalar == F::ZERO {
				continue;
			}

			// otherwise, we permute both our matrix and our permutation
			// to get the stable_index to be in the appropriate row.
			perm.swap(r, stable_index);
			mat.permute(r, stable_index);

			// permute L matrix, but only entries under the diagonal
			if c>0{
				for c2 in 0..c {
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
				let scalar = mat[(r2,c)];
				mat.replace(r, r2, F::ZERO-scalar);
				elem_prod[(r2, r)] = scalar;
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

	#[test]
	fn test_lu(){
		let mat : Mat<f32> = Mat::new(3,4,
			vec![
				0.0, 1.0, 2.0,  3.0,
				4.0, 5.0, 6.0,  7.0,
				8.0, 9.0, 10.0, 11.0
			]
		);
	}
}