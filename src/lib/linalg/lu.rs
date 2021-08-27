use super::super::num::{Field, StabilityCmp, EpsilonEquality};
use super::mat::Mat;
use super::util::get_max_index;


pub struct RREFDecomposition<F: Field> {
	pub p : Mat<F>,
	pub rref : Mat<F>
}

pub struct LUDecomposition<F: Field> {
	pub p : Mat<F>,
	pub l : Mat<F>,
	pub u : Mat<F>
}

struct LUDecompositionInternal<F: Field> {
	pub p: Vec<usize>,
	pub l: Mat<F>,
	pub u: Mat<F>
}

impl<F: Field> LUDecompositionInternal<F> {

	// assumes src < target
	fn permute(&mut self, src: usize, target: usize) {
		self.p.swap(src, target);
		self.u.permute(src, target);
		self.l.permute_under_diagonal(src, target);
	}

	fn scale_row(&mut self, r: usize, scalar: F) {
		self.l[(r,r)] = scalar;
		self.u.scale(r, F::ONE/scalar);
	}

	fn replace_col_under_row(&mut self, r: usize, c: usize) {
		for r2 in (r+1)..self.u.rows() {
			let scalar2 = self.u[(r2,c)];
			self.u.replace(r, r2, F::ZERO-scalar2);
			self.l[(r2, r)] = scalar2;
		}
	}
}

impl<F: Field> Mat<F> {
	// assumes src < target
	fn permute_under_diagonal(&mut self, src: usize, target: usize) {
		for c in 0..src {
			let temp = self[(src, c)];
			self[(src, c)] = self[(target, c)];
			self[(target, c)] = temp;
		}
	}
}

impl<F: Field+StabilityCmp+EpsilonEquality> Mat<F> {

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

	/// Returns LU Decomposition of given matrix. This method does not make
	/// any assumptions on the dimension of the matrix (i.e. we allow
	/// rectangular matrices), though it does require the field to implement 
	/// StabilityCmp in order to do Gaussian elimination under the hood.
	/// 
	/// The LU decomposition satisfies more requirements than that typically
	/// seen in Fortran implementations. The `u` matrix is actually in row
	/// echelon form, which will be convenient for calculating rank as well as
	/// rref (which will be used to compute a basis for the kernel)
	/// 
	/// `l` is product of inverses of elementary matrices applied to `mat`,
	/// `p` is permutation, such that `(p * self).epsilon_equals(l * u)`
	/// or equivalently `self.epsilon_equals(p.transpose() * l * u)`.
	pub fn lu(&self) -> LUDecomposition<F> {
		let m = self.cols();
		let n = self.rows();

		// the "in progress" LU decomposition that will eventually become what
		// we want.
		let mut lu = LUDecompositionInternal {
			p: (0..n).collect(),
			l: Mat::identity(n),
			u: self.clone()
		};

		// row that hasn't been cleared yet
		let mut r = 0;

		for c in 0..m {
			if r >= n {
				break;
			}

			if let Some((stable_index, scalar)) = lu.u.get_stable_column_index_under_row(c, r) {
				lu.permute(r, stable_index);
				lu.scale_row(r, scalar);
				lu.replace_col_under_row(r, c);
				
				r += 1;
			}
		}
		LUDecomposition{
			p: Mat::permutation_from_vec(lu.p),
			l: lu.l,
			u: lu.u
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