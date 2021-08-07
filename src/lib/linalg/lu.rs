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
		let mut perm : Vec<usize> = (0..n-1).collect();
		let mut elem_prod : Mat<F> = Mat::identity(n);

		for c in 0..mat.cols()-1 {
			let current_col = mat.get_col(c);

			// get index of most stable (usually largest) entry in column
			let stable_index = get_max_index(
				current_col.entries()[c..].iter(),
				|x,y| x.stability_cmp(y)
			)+c;

			// if "most stable" index is zero, then no replacement is needed
			// in this column, so we skip to the next.
			if mat[(stable_index, c)]==F::ZERO {
				continue;
			} 

			// otherwise, we permute both our matrix and our permutation
			// to get the stable_index to be in the appropriate row.
			mat.permute(c, stable_index);
			mat.scale(c, F::ONE/mat[(c,c)]);
			//mat = &Mat::scale(n, c, F::ONE/ *current_col.get_unchecked(stable_index,0)) * &(&temp_perm * &mat);
			perm.swap(c, stable_index);

			// SECOND: do replacement
			for r in (c+1)..mat.rows() {
				let scalar = mat[(r,c)];
				mat.replace(c, r, F::ZERO-scalar);
				elem_prod[(c, r)] = scalar;
			}
			
		}
		LUDecomposition{
			p: Mat::permutation_from_vec(perm),
			l: mat,
			u: elem_prod
		}
	}
}