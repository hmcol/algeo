use super::super::num::{Field, StabilityCmp};
use super::mat::Mat;
use super::util::get_max_index;
use std::collections::HashSet;

pub struct LUDecomposition<F: Field> {
	p : Mat<F>,
	l : Mat<F>,
	u : Mat<F>
}

impl<F: Field+StabilityCmp> Mat<F> {

	/// NOT IMPLEMENTED YET. THIS IS REMNANTS FROM ROW ECHELON FORM METHOD.
	pub fn lu(&self) -> LUDecomposition<F> {
		// Matrix that is eventually made into row echelon form. It
		// will become upper triangular.
		let mut mat = self.clone();
		let n = mat.rows();

		// elem_prod is product of inverses of elementary matrices applied
		// to `mat`, perm is permutation matrix, such that
		// `perm * self=elem_prod * mat`
		// or equivalently
		// `self = perm^T * elem_prod * mat`
		let mut perm : Mat<F> = Mat::identity(n);
		let mut elem_prod : Mat<F> = Mat::identity(n);

		for c in 0..mat.cols()-1 {
			let current_col = mat.get_col(c);

			// If not all entries below are zero
			if !current_col.entries()[c..].iter().all(|x| *x==F::ZERO) {
				// FIRST: put best row at the top

				// get index of most stable (usually largest) entry in column
				let stable_index = get_max_index(
					current_col.entries()[c..].iter(),
					|x,y| x.stability_cmp(y)
				)+c;

				let temp_perm = Mat::permutation(n, c, stable_index);
				mat = &Mat::scale(n, c, F::ONE/ *current_col.get_unchecked(stable_index,0)) * &(&temp_perm * &mat);
				perm = &temp_perm * &perm;

				// SECOND: do replacement
				for r in (c+1)..mat.rows() {
					mat = &Mat::replacement(n, c, r, F::ZERO - *mat.get_unchecked(r,c)) * &mat;
				}
			}
		}
		LUDecomposition{
			p: perm,
			l: mat,
			u: elem_prod
		}
	}
}