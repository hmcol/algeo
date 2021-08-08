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
				mat.get_col_iter(c).skip(std::cmp::max(0, (r as i32) -1) as usize),
				|x,y| x.stability_cmp(y)
			)+r;

			// if "most stable" index is zero, then no replacement is needed
			// in this column, so we skip to the next.
			if mat[(stable_index, c)]==F::ZERO {
				continue;
			}

			// otherwise, we permute both our matrix and our permutation
			// to get the stable_index to be in the appropriate row.
			mat.permute(r, stable_index);
			mat.scale(r, F::ONE/mat[(c,c)]);
			perm.swap(r, stable_index);

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