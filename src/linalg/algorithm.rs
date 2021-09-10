use std::{collections::HashSet, ops::IndexMut};

use super::mat::Mat;
use super::util::get_max_index;
use crate::core::num::{EpsilonEquality, Field, StabilityCmp};

// utility matrix methods ------------------------------------------------------

impl<F: Field + StabilityCmp + EpsilonEquality> Mat<F> {
    /// looks at column `c`, and returns the row under (or equal) to `r` such
    /// that the entry `self[(r,c)]` is maximal under the ordering given by
    /// the `StabilityCmp` trait. In practice, this method is used to a good
    /// pivot. For example, we want to avoid dividing by `f32`s of small
    /// magnitude, like `0.00001`, since this will cause more floating point
    /// error.
    pub fn get_stable_column_index_under_row(&self, c: usize, r: usize) -> Option<(usize, F)> {
        let stable_index =
            get_max_index(self.get_col_iter(c).skip(r), |x, y| x.stability_cmp(y)) + r;

        let scalar = self[(stable_index, c)];
        if scalar.epsilon_equals(&F::ZERO) {
            None
        } else {
            Some((stable_index, scalar))
        }
    }
}

impl<F: Field + EpsilonEquality> Mat<F> {
    /// Looks at the row `r` and returns the first column (left-to-right) such
    /// that `self[(r,c)]` is nonzero (up to `epsilon_equals`)
    fn index_of_first_nonzero_entry(&self, r: usize) -> Option<usize> {
        self.get_row(r)
            .entries()
            .iter()
            .position(|&x| !x.epsilon_equals(&F::ZERO))
    }

    /// Counts the number of zero rows (up to `epsilon_equals`).
    pub fn zero_rows(&self) -> usize {
        let mut zero_rows = 0;

        for r in 0..self.rows() {
            if let None = self.index_of_first_nonzero_entry(r) {
                zero_rows += 1;
            }
        }

        zero_rows
    }
}

impl<F: Field> Mat<F> {
    // assumes src < target
    fn permute_under_diagonal(&mut self, src: usize, target: usize) {
        for c in 0..src {
            self.swap((src, c), (target, c));
            // let temp = self[(src, c)];
            // self[(src, c)] = self[(target, c)];
            // self[(target, c)] = temp;
        }
    }
}

// LU Decomposition ------------------------------------------------------------

pub struct LUDecomposition<F: Field> {
    pub p: Mat<F>,
    pub l: Mat<F>,
    pub u: Mat<F>,
}

struct LUDecompositionInternal<F: Field> {
    pub p: Vec<usize>,
    pub l: Mat<F>,
    pub u: Mat<F>,
}

impl<F: Field> LUDecompositionInternal<F> {
    // assumes src < target
    fn permute(&mut self, src: usize, target: usize) {
        self.p.swap(src, target);
        self.u.permute_rows(src, target);
        self.l.permute_under_diagonal(src, target);
    }

    fn scale_row(&mut self, r: usize, scalar: F) {
        self.l[(r, r)] = scalar;
        self.u.scale_row(r, F::ONE / scalar);
    }

    fn replace_col_under_row(&mut self, r: usize, c: usize) {
        for r2 in (r + 1)..self.u.rows() {
            let scalar2 = self.u[(r2, c)];
            self.u.replace_row(r, r2, F::ZERO - scalar2);
            self.l[(r2, r)] = scalar2;
        }
    }
}

impl<F: Field + StabilityCmp + EpsilonEquality> Mat<F> {
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
            u: self.clone(),
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
        LUDecomposition {
            p: Mat::permutation_from_vec(lu.p),
            l: lu.l,
            u: lu.u,
        }
    }
}

// Row Equivalent Forms --------------------------------------------------------

/// Intended to be used as pair of invertible matrix `p` and row equivalent
/// matrix `b` such that `p*a=b`. In the particular case that `b` is the
/// identity, `p` is the inverse of `a`.
pub trait RowEquivalentForm<F: Field> {
    fn p(&mut self) -> &mut Mat<F>;
    fn b(&mut self) -> &mut Mat<F>;

    fn permute(&mut self, r1: usize, r2: usize) {
        self.p().permute_rows(r1, r2);
        self.b().permute_rows(r1, r2);
    }

    fn scale(&mut self, r: usize, scalar: F) {
        self.p().scale_row(r, scalar);
        self.b().scale_row(r, scalar);
    }

    /// Performs replacement operation to clear all entries in column under the
    /// given row. Assumes `self[(r,c)] == F::ONE`.
    fn replace_col_under_row(&mut self, r: usize, c: usize) {
        for r2 in (r + 1)..self.b().rows() {
            let scalar2 = self.b()[(r2, c)];
            self.p().replace_row(r, r2, F::ZERO - scalar2);
            self.b().replace_row(r, r2, F::ZERO - scalar2);
        }
    }

    /// Performs replacement operation to clear all entries in column above the
    /// given row. Assumes `self[(r,c)] == F::ONE`.
    fn replace_col_above_row(&mut self, r: usize, c: usize) {
        for r2 in 0..r {
            let scalar2 = self.b()[(r2, c)];
            self.p().replace_row(r, r2, F::ZERO - scalar2);
            self.b().replace_row(r, r2, F::ZERO - scalar2);
        }
    }
}

/// Struct containing original matrix `original`, row echelon form `b`, and
/// product of row operations `p`.
pub struct RowEchelonForm<'a, F: Field> {
    pub p: Mat<F>,
    pub b: Mat<F>,
    pub original: &'a Mat<F>,
}

impl<'a, F: Field> RowEquivalentForm<F> for RowEchelonForm<'a, F> {
    fn p(&mut self) -> &mut Mat<F> {
        &mut self.p
    }

    fn b(&mut self) -> &mut Mat<F> {
        &mut self.b
    }
}

/// Struct containing original matrix `original`, reduced row echelon form `b`,
/// and product of row operations `p`.
pub struct ReducedRowEchelonForm<'a, F: Field> {
    pub p: Mat<F>,
    pub rref: Mat<F>,
    pub original: &'a Mat<F>,
}

impl<'a, F: Field> RowEquivalentForm<F> for ReducedRowEchelonForm<'a, F> {
    fn p(&mut self) -> &mut Mat<F> {
        &mut self.p
    }

    fn b(&mut self) -> &mut Mat<F> {
        &mut self.rref
    }
}

impl<F: Field + StabilityCmp + EpsilonEquality> Mat<F> {
    /// Computes row echelon form
    pub fn row_echelon(&self) -> RowEchelonForm<F> {
        let n = self.rows();
        let m = self.cols();

        let mut form = RowEchelonForm {
            b: self.clone(),
            p: Mat::identity(n),
            original: &self,
        };

        // row that hasn't been cleared yet
        let mut r = 0;

        for c in 0..m {
            if r >= n {
                break;
            }

            if let Some((stable_index, scalar)) = form.b.get_stable_column_index_under_row(c, r) {
                form.permute(r, stable_index);
                form.scale(r, F::ONE / scalar);
                form.replace_col_under_row(r, c);

                r += 1;
            }
        }

        form
    }
}

impl<'a, F: Field + StabilityCmp + EpsilonEquality> RowEchelonForm<'a, F> {
    /// Computes the nullity (dimension of the kernel)
    pub fn nullity(&self) -> usize {
        self.b.cols() - self.rank()
    }

    /// Computes the rank (dimension of the range)
    pub fn rank(&self) -> usize {
        self.b.rows() - self.b.zero_rows()
    }

    /// Computes reduced row echelon form
    pub fn to_rref(self) -> ReducedRowEchelonForm<'a, F> {
        let n = self.b.rows();

        let mut form = ReducedRowEchelonForm {
            p: self.p,
            rref: self.b,
            original: self.original,
        };

        for r in (1..n).rev() {
            if let Some(c) = form.rref.index_of_first_nonzero_entry(r) {
                form.replace_col_above_row(r, c);
            }
        }

        form
    }
}

impl<'a, F: Field + StabilityCmp + EpsilonEquality> ReducedRowEchelonForm<'a, F> {
    /// Computes the nullity (dimension of the kernel)
    pub fn nullity(&self) -> usize {
        self.rref.cols() - self.rank()
    }

    /// Computes the rank (dimension of the range)
    pub fn rank(&self) -> usize {
        self.rref.rows() - self.rref.zero_rows()
    }

    /// Computes basis for the kernel
    pub fn compute_kernel(&self) -> Vec<Mat<F>> {
        let n = self.rref.cols();
        let bound_variables: HashSet<usize> = (0..self.rref.rows())
            .filter_map(|r| self.rref.index_of_first_nonzero_entry(r))
            .collect();

        let mut basis: Vec<Mat<F>> = vec![];

        for c in (0..self.rref.cols()).filter(|c| !bound_variables.contains(&c)) {
            // if free variable
            let mut vector: Mat<F> = Mat::e(c, n);
            for r in 0..self.rref.rows() {
                if let Some(c2) = self.rref.index_of_first_nonzero_entry(r) {
                    if c2 < c {
                        *vector.get_mut_unchecked(c2, 0) = F::ZERO - *self.rref.get_unchecked(r, c);
                    }
                }
            }
            basis.push(vector);
        }

        basis
    }

    /// Computes basis for the range. UNTESTED.
    pub fn compute_range(&self) -> Vec<Mat<F>> {
        let bound_variables: HashSet<usize> = (0..self.rref.rows())
            .filter_map(|r| self.rref.index_of_first_nonzero_entry(r))
            .collect();

        bound_variables
            .into_iter()
            .map(|c| self.original.get_col(c))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use crate::core::num::Field;
    use crate::linalg::util::mat_iterator;

    use super::Mat;
    use crate::core::num::EpsilonEquality;

    #[test]
    fn test_lu() {
        let n: usize = 3;
        let m: usize = 4;
        let values = vec![0.0, 1.0, 2.0];

        // iterate through all matrices with only entries in `values`
        for mat in mat_iterator(n, m, &values) {
            let lu = mat.lu();
            let prod = &lu.p.transpose() * &(&lu.l * &lu.u);

            assert!(
                lu.l.is_lower_triangular(),
                "L={} is not lower triangular, original A={}",
                lu.l,
                mat
            );
            assert!(
                lu.u.is_upper_triangular(),
                "U={} is not upper triangular, original A={}",
                lu.u,
                mat
            );
            assert!(
                mat.epsilon_equals(&prod),
                "failed on {} == {}, epsilon error is {}",
                mat,
                prod,
                (&mat + &(&prod * (-1.0))).frobenius_norm()
            );
        }
    }

    #[test]
    fn test_row_echelon() {
        let n: usize = 3;
        let m: usize = 4;
        let values = vec![0.0, 1.0, 2.0];

        // iterate through all matrices with only entries in `values`
        for mat in mat_iterator(n, m, &values) {
            let row_ech = mat.row_echelon();
            let lu = mat.lu();

            assert!(
                lu.u.epsilon_equals(&row_ech.b),
                "failed on {} == {}, epsilon error is {}, does not match LU decomposition",
                lu.u,
                row_ech.b,
                (&lu.u + &(&row_ech.b * (-1.0))).frobenius_norm()
            );

            let prod = &row_ech.p * &mat;
            assert!(
                prod.epsilon_equals(&row_ech.b),
                "failed on {} == {}, epsilon error is {}, does not satisfy `p*a=b`",
                prod,
                row_ech.b,
                (&prod + &(&row_ech.b * (-1.0))).frobenius_norm()
            );
        }
    }

    impl<F: Field + EpsilonEquality> Mat<F> {
        fn is_col_cleared(&self, c: usize, pivot_r: usize) -> bool {
            for r in 0..self.rows() {
                if r != pivot_r && !self[(r, c)].epsilon_equals(&F::ZERO) {
                    return false;
                }
            }

            true
        }

        fn is_rref(&self) -> bool {
            let mut prev_col: i32 = -1;

            for r in 0..self.rows() {
                if let Some(c) = self.index_of_first_nonzero_entry(r) {
                    if !self.is_col_cleared(c, r) {
                        return false;
                    }

                    if prev_col >= c as i32 {
                        return false;
                    }
                    prev_col = c as i32;
                } else {
                    prev_col = self.cols() as i32;
                }
            }

            true
        }
    }

    #[test]
    fn test_rref() {
        let n: usize = 3;
        let m: usize = 4;
        let values = vec![0.0, 1.0, 2.0];

        // iterate through all matrices with only entries in `values`
        for mat in mat_iterator(n, m, &values) {
            let rref = mat.row_echelon().to_rref();

            let prod = &rref.p * &mat;
            assert!(
                prod.epsilon_equals(&rref.rref),
                "failed on {} == {}, epsilon error is {}, does not satisfy `p*a=b`",
                prod,
                rref.rref,
                (&prod + &(&rref.rref * (-1.0))).frobenius_norm()
            );

            assert!(rref.rref.is_rref(), "is not rref");
        }
    }

    #[test]
    fn test_compute_kernel() {
        let n: usize = 3;
        let m: usize = 4;
        let values: Vec<f64> = vec![0.0, 1.0, 2.0];

        // iterate through all matrices with only entries in `values`
        for mat in mat_iterator(n, m, &values) {
            let row_ech = mat.row_echelon();
            let nullity = row_ech.nullity();

            let rref = row_ech.to_rref();
            let kernel_basis = rref.compute_kernel();

            // test that vectors are in the kernel
            for v in kernel_basis.iter() {
                let prod = &mat * v;
                assert!(
                    prod.epsilon_equals(&Mat::zero(prod.rows(), prod.cols())),
                    "failed on {} not in kernel, norm is {}",
                    prod,
                    prod.frobenius_norm()
                );
            }

            // check kernel vectors are linearly independent
            let kernel_mat = Mat::from_row_vectors(&kernel_basis);
            let ker_row_ech = kernel_mat.row_echelon();
            assert_eq!(ker_row_ech.rank(), nullity);
        }
    }
}
