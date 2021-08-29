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

	/// TODO: Needs testing, it very likely does not work
	/// Input: generators for the cone
	/// Output: Cone<F>, which computes from the generators a list of facets and
	/// a basis for the orthogonal complement.
	pub fn new(generators: Vec<Mat<F>>) -> Cone<F> {
		let ambient_dim = generators[0].rows();

		let generator_mat = Mat::from_row_vectors(&generators);
		let rref = generator_mat.row_echelon().to_rref();
		let dim = rref.rank();

		let orth_comp_basis = rref.compute_kernel();

		let mut facets: Vec<Face<F>> = vec![];

		// compute facets
		'outer: for subset in generators.iter().combinations(dim-1) {
			// this is very inefficient to copy this all the time
			let subset_owned : Vec<Mat<F>> = subset.iter().map(|&m| m.clone()).collect();
			let subset_mat = Mat::from_row_vectors(&subset_owned);
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

				// re-orient the normal in the direction of the first generator
				let is_above = match ordering {
					std::cmp::Ordering::Less => false,
					_ => true
				};
				if !is_above {
					normal.scale(F::ZERO-F::ONE);
				}

				// check that the remaining generators "are above" normal
				for g in generators.iter() {
					if let Some(ordering) = g.dot(&normal).partial_cmp(&F::ZERO) {
						// TODO add clause checking for when dot product epsilon_equals 0
						// in that case, add that generator to the subset (or make a new collection)
						// of extra vectors. This will give a way for combining "facets" that are
						// identical.
						match ordering {
							std::cmp::Ordering::Less => continue 'outer,
							_ => {}
						};
					} else {
						continue 'outer; 
					}
				}

				// create facet to add to list
				let facet = Face {
					normal,
					generators: subset_owned
				};

				facets.push(facet);
			}
		}

		Cone {
			generators,
			facets,
			orth_comp_basis,
			ambient_dim,
			dim
		}
	}
}