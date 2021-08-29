use itertools::Itertools;

use crate::lib::num::{EpsilonEquality, Field, StabilityCmp};


use super::mat::Mat;


pub struct Cone<F: Field> {
	generators: Vec<Mat<F>>,
	facets: Vec<Face<F>>,
	orth_comp_basis: Vec<Mat<F>>,

	ambient_dim: usize,
	dim: usize,
}

pub struct Face<F: Field> {
	normal: Mat<F>,
	generators: Vec<Mat<F>>
}

impl<F: Field + StabilityCmp + EpsilonEquality + PartialOrd> Cone<F> {

	pub fn new(generators: Vec<Mat<F>>) -> Cone<F> {
		let ambient_dim = generators[0].rows();

		let generator_mat = Mat::from_row_vectors(&generators);
		let rref = generator_mat.row_echelon().to_rref();
		let dim = rref.rank();

		let orth_comp_basis = rref.compute_kernel();

		// compute facets
		for subset in generators.iter().combinations(dim-1) {
			// this is very inefficient to copy this all the time
			let owned_subset : Vec<Mat<F>> = subset.iter().map(|&m| m.clone()).collect();
			let subset_mat = Mat::from_row_vectors(&owned_subset);
			let subset_rref = subset_mat.row_echelon().to_rref();

			// (1) check linearly independent
			if subset_rref.rank() != dim-1 {
				continue;
			}

			// (2) get a normal (there are potentially multiple in dim != ambient_dim)
			let mut subset_kernel = subset_rref.compute_kernel();
			let mut normal = subset_kernel[0].clone();

			// (3) check that all generators are on one side
			if let Some(ordering) = generators[0].dot(&normal).partial_cmp(&F::ZERO) {
				let is_above = match ordering {
					std::cmp::Ordering::Less => false,
					_ => true
				};
				if !is_above {
					normal.scale_row(1, F::ZERO-F::ONE);
				}

				for g in generators.iter() {
					
				}
			}
		}

		Cone {
			generators,
			facets: todo!(),
			orth_comp_basis,
			ambient_dim,
			dim
		}
	}
}