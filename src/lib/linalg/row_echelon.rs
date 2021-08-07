

use super::super::num::{Field, StabilityCmp};
use super::mat::Mat;
use super::util::get_max_index;
use std::collections::HashSet;

impl<F: Field+StabilityCmp> Mat<F> {

	pub fn row_echelon(&self) -> (Mat<F>, Mat<F>) {
		let mut temp = self.clone();
		let n = temp.rows();
		let mut permutation : Mat<F> = Mat::identity(n);

		for c in 0..temp.cols()-1 {
			let current_col = temp.get_col(c);

			// if not all zero
			if !current_col.entries()[c..].iter().all(|x| *x==F::ZERO) {
				// FIRST: put best row at the top

				// get index of most stable (usually largest) entry in column
				let stable_index = get_max_index(
					current_col.entries()[c..].iter(),
					|x,y| x.stability_cmp(y)
				)+c;

				let temp_perm = Mat::permutation(n, c, stable_index);
				temp = &Mat::scale_mat(n, c, F::ONE/ *current_col.get_unchecked(stable_index,0)) * &(&temp_perm * &temp);
				permutation = &temp_perm * &permutation;

				// SECOND: do replacement
				for r in (c+1)..temp.rows() {
					temp = &Mat::replacement(n, c, r, F::ZERO - *temp.get_unchecked(r,c)) * &temp;
				}
			}
		}
		(temp, permutation)
	}

	pub fn index_of_first_nonzero_entry(&self, r: usize) -> Option<usize> {
		self.get_row(r).entries().iter().position(|x| *x != F::ZERO)
	}

	pub fn reduced_row_echelon(&self) -> (Mat<F>, Mat<F>) {
		let (row_echelon, permutation) = self.row_echelon();
		let mut temp = row_echelon;
		let n = temp.rows();

		for r in (1..n).rev() {
			let c_option = temp.index_of_first_nonzero_entry(r);
			if let Some(c) = c_option {
				for r2 in 0..r {
					temp = &Mat::replacement(n, r, r2, F::ZERO-*temp.get_unchecked(r2,c)) * &temp;
				}
			}
		}

		(temp, permutation)
	}

	// computes a basis for the kernel of matrix
	pub fn compute_kernel(&self) -> Vec<Mat<F>> {
		let (rref, permutation) = self.reduced_row_echelon();
		let n = rref.cols();
		let bound_variables : HashSet<usize> = (0..rref.rows()).filter_map(|r| rref.index_of_first_nonzero_entry(r)).collect();

		let mut basis : Vec<Mat<F>> = vec![];

		for c in (0..rref.cols()).filter(|c| !bound_variables.contains(&c)) {
			// if free variable
			let mut vector : Mat<F> = Mat::e(c, n);
			for r in 0..rref.rows() {
				if let Some(c2) = rref.index_of_first_nonzero_entry(r) {
					if c2 < c {
						*vector.get_mut_unchecked(c2, 0) = F::ZERO-*rref.get_unchecked(r,c);
					}
				}
			}
			basis.push(vector);
		}

		basis
	}
}