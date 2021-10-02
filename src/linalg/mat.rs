use itertools::Itertools;
use xops::binop;

use super::util::get_box_iter;
use crate::core::num::{EpsilonEquality, Field};
use std::ops::{Add, Index, IndexMut, Mul, MulAssign};
use std::fmt;

#[derive(Debug, PartialEq, Clone)]
pub struct Mat<F: Field> {
    entries: Vec<F>,
    rows: usize,
    cols: usize,
}

impl<F: Field> Mat<F> {
    pub fn new(rows: usize, cols: usize, entries: Vec<F>) -> Mat<F> {
        Mat {
            entries,
            rows,
            cols,
        }
    }

    pub fn swap(&mut self, (r1, c1): (usize, usize), (r2, c2): (usize, usize)) {
        self.entries.swap(r1 * self.cols + c1, r2 * self.cols + c2);
    }

    // Mat<F> constructors ----------------------------------------------------

    /// identity matrix with of dimension `n * n`
    pub fn identity(n: usize) -> Mat<F> {
        Mat::new(
            n,
            n,
            get_box_iter(n, n)
                .map(|(r, c)| if r == c { F::ONE } else { F::ZERO })
                .collect(),
        )
    }

    pub fn from_col_vectors(vs: &[Mat<F>]) -> Mat<F> {
        Mat::new(
            vs[0].rows(),
            vs.len(),
            get_box_iter(vs[0].rows(), vs.len())
                .map(|(r, c)| vs[c][(r, 0)])
                .collect(),
        )
    }

    pub fn from_row_vectors(vs: &[Mat<F>]) -> Mat<F> {
        Mat::new(
            vs.len(),
            vs[0].rows(),
            get_box_iter(vs.len(), vs[0].rows())
                .map(|(r, c)| vs[r][(c, 0)])
                .collect(),
        )
    }

    /// `k`th standard basis vector of dimension `n`
    pub fn e(k: usize, n: usize) -> Mat<F> {
        Mat::new(
            n,
            1,
            (0..n)
                .map(|r| if r == k { F::ONE } else { F::ZERO })
                .collect(),
        )
    }

    /// zero matrix of dimension `rows * cols`
    pub fn zero(rows: usize, cols: usize) -> Mat<F> {
        Mat::new(
            rows,
            cols,
            get_box_iter(rows, cols).map(|(_r, _c)| F::ZERO).collect(),
        )
    }

    /// elementary matrix of dimension `rows*rows` that scales the `i`th row
    pub fn scale_mat(rows: usize, i: usize, scale: F) -> Mat<F> {
        Mat::new(
            rows,
            rows,
            get_box_iter(rows, rows)
                .map(|(r, c)| {
                    if r == c && r == i {
                        scale
                    } else if r == c {
                        F::ONE
                    } else {
                        F::ZERO
                    }
                })
                .collect(),
        )
    }

    /// elementary matrix of dimension `rows*rows` that permutes the `source`
    /// and `target` rows.
    pub fn permutation(rows: usize, source: usize, target: usize) -> Mat<F> {
        Mat::new(
            rows,
            rows,
            get_box_iter(rows, rows)
                .map(|(r, c)| {
                    if (r, c) == (source, target)
                        || (c, r) == (source, target)
                        || (r == c && r != source && c != target)
                    {
                        F::ONE
                    } else {
                        F::ZERO
                    }
                })
                .collect(),
        )
    }

    /// input - permutation : Vec<F> of length n with entries `0..(n-1)` in
    /// some order. The vec can be thought of as a bijection of
    /// {0,1,...,n-1} -> {0,1,...,n-1}
    /// and `permutation_from_vec` returns a matrix that does this permutation
    /// to the rows.
    pub fn permutation_from_vec(permutation: Vec<usize>) -> Mat<F> {
        let n = permutation.len();
        Mat::new(
            n,
            n,
            get_box_iter(n, n)
                .map(|(r, c)| if permutation[r] == c { F::ONE } else { F::ZERO })
                .collect(),
        )
    }

    /// elementary matrix of dimension `rows*rows` adds the `source` row
    /// scaled by `scalar` to the `target` row.
    pub fn replacement(rows: usize, source: usize, target: usize, scalar: F) -> Mat<F> {
        Mat::new(
            rows,
            rows,
            get_box_iter(rows, rows)
                .map(|(r, c)| {
                    if c == source && r == target {
                        scalar
                    } else if r == c {
                        F::ONE
                    } else {
                        F::ZERO
                    }
                })
                .collect(),
        )
    }

    /// returns transposition of the matrix.
    pub fn transpose(&self) -> Mat<F> {
        Mat::new(
            self.cols,
            self.rows,
            get_box_iter(self.cols, self.rows)
                .map(|(c, r)| self.get_unchecked(r, c))
                .cloned()
                .collect(),
        )
    }

    // Mat<F> mutators ---------------------------------------------------------

    pub fn scale_row(&mut self, i: usize, scalar: F) {
        for c in 0..self.cols {
            self[(i, c)] = scalar * self[(i, c)];
        }
    }

    pub fn permute_rows(&mut self, source: usize, target: usize) {
        for c in 0..self.cols {
            let temp = self[(source, c)];
            self[(source, c)] = self[(target, c)];
            self[(target, c)] = temp;
        }
    }

    pub fn replace_row(&mut self, source: usize, target: usize, scalar: F) {
        for c in 0..self.cols {
            self[(target, c)] = self[(target, c)] + self[(source, c)] * scalar;
        }
    }

    // TODO Refactor into MulAssign
    pub fn scale(&mut self, scalar: F) {
        self.entries.iter_mut().for_each(|x| *x = *x * scalar);
    }

    // Mat<F> getters (possibly mutable) ---------------------------------------

    pub fn get_unchecked(&self, r: usize, c: usize) -> &F {
        &self.entries[r * self.cols + c]
    }

    pub fn entries(&self) -> &[F] {
        &self.entries
    }

    pub fn get_mut_unchecked(&mut self, r: usize, c: usize) -> &mut F {
        &mut self.entries[r * self.cols + c]
    }

    pub fn rows(&self) -> usize {
        self.rows
    }
    pub fn cols(&self) -> usize {
        self.cols
    }

    pub fn get_row_iter(&self, r: usize) -> impl Iterator<Item = &F> {
        self.entries[r * self.cols..(r + 1) * self.cols].iter()
    }

    pub fn get_row_iter_mut(&mut self, r: usize) -> impl Iterator<Item = &mut F> {
        self.entries[r * self.cols..(r + 1) * self.cols].iter_mut()
    }

    /// transposed to be column vector
    pub fn get_row(&self, r: usize) -> Mat<F> {
        Mat::from(self.get_row_iter(r).cloned())
    }

    pub fn get_col_iter(&self, c: usize) -> impl Iterator<Item = &F> {
        (0..self.rows).map(move |r| self.get_unchecked(r, c))
    }

    pub fn get_col_iter_mut(&mut self, c: usize) -> impl Iterator<Item = &mut F> {
        // magic code that I got from discord.
        // the more simple (0..self.rows).map(move |r| entries.get_mut(r,c))
        // has weird lifetime issues. This works though, so sicko mode
        self.entries
            .chunks_exact_mut(self.rows)
            .map(move |row| &mut row[c])
    }

    pub fn get_col(&self, c: usize) -> Mat<F> {
        Mat::from(self.get_col_iter(c).cloned())
    }

    pub fn frobenius_norm(&self) -> F {
        self.entries.iter().fold(F::ZERO, |accum, &x| accum + x * x)
    }

    // Operations --------------------------------------------------------------

    pub fn dot(&self, other: &Mat<F>) -> F {
        self.entries
            .iter()
            .zip(other.entries())
            .fold(F::ZERO, |accum, (&x, &y)| accum + x * y)
    }

    // Util --------------------------------------------------------------------

    pub fn index_iter(&self) -> impl Iterator<Item = (usize, usize)> {
        (0..self.rows()).cartesian_product(0..self.cols())
    }

    // Debug/Testing methods ---------------------------------------------------

    /// TODO: decide whether `is_upper_triangular` should require F to implement
    /// `EpsilonEquality`. In practice, I'm only using `is_upper_triangular` to
    /// verify the correctness of the LU decomposition, and in that case the
    /// replacement operations are exact (so that I wouldn't have to loosen
    /// equality to epsilon equality).
    pub fn is_upper_triangular(&self) -> bool {
        for r in 0..self.rows {
            for c in 0..(std::cmp::min(r, self.cols)) {
                if self[(r, c)] != F::ZERO {
                    return false;
                }
            }
        }

        true
    }

    pub fn is_lower_triangular(&self) -> bool {
        for r in 0..self.rows {
            for c in (std::cmp::min(r, self.cols) + 1)..self.cols {
                if self[(r, c)] != F::ZERO {
                    return false;
                }
            }
        }

        true
    }
}

impl<F: Field, I: Iterator<Item = F>> From<I> for Mat<F> {
    fn from(iter: I) -> Self {
        let entries: Vec<F> = iter.collect();
        let rows = entries.len();
        Mat {
            entries,
            rows,
            cols: 1,
        }
    }
}

impl<F: Field> Index<(usize, usize)> for Mat<F> {
    type Output = F;

    fn index(&self, (r, c): (usize, usize)) -> &F {
        self.get_unchecked(r, c)
    }
}

impl<F: Field> IndexMut<(usize, usize)> for Mat<F> {
    fn index_mut(&mut self, (r, c): (usize, usize)) -> &mut F {
        self.get_mut_unchecked(r, c)
    }
}

impl<F: Field> fmt::Display for Mat<F> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "\n")?;
        for r in 0..self.rows {
            f.debug_list()
                .entries(&self.entries[r * self.cols..r * self.cols + self.cols])
                .finish()?;
            if r != self.rows - 1 {
                write!(f, "\n")?;
            }
        }
        Ok(())
    }
}

#[binop(derefs)]
impl<F: Field> Add for &Mat<F> {
    type Output = Mat<F>;

    fn add(self, other: &Mat<F>) -> Mat<F> {
        if !(self.rows == other.rows && self.cols == other.cols) {
            panic!("tried to add matrix with incompatible dimensions.")
        }
        return Mat::new(
            self.rows,
            self.cols,
            self.entries
                .iter()
                .zip(other.entries.iter())
                .map(|(&x, &y)| x + y)
                .collect::<Vec<F>>(),
        );
    }
}

// scalar multiplication
#[binop(derefs)]
impl<F: Field> Mul<F> for &Mat<F> {
    type Output = Mat<F>;

    fn mul(self, scalar: F) -> Mat<F> {
        return Mat::new(
            self.rows,
            self.cols,
            self.entries.iter().map(|&x| x * scalar).collect(),
        );
    }
}

#[binop(derefs)]
impl<F: Field> Mul for &Mat<F> {
    type Output = Mat<F>;

    fn mul(self, other: &Mat<F>) -> Mat<F> {
        if self.cols != other.rows {
            panic!("tried to multiply matrices with incompatible dimensions.")
        }
        Mat::new(
            self.rows,
            other.cols,
            get_box_iter(self.rows, other.cols)
                .map(|(r, c)| {
                    self.get_row_iter(r)
                        .zip(other.get_col_iter(c))
                        .fold(F::ZERO, |accum: F, (&x, &y): (&F, &F)| accum + x * y)
                })
                .collect(),
        )
    }
}

impl<F: Field> MulAssign<F> for Mat<F> {
    fn mul_assign(&mut self, scalar: F) {
		self.entries.iter_mut().for_each(|x| *x *= scalar);
    }
}

impl<F: Field + EpsilonEquality> EpsilonEquality for Mat<F> {
    fn epsilon_equals(&self, other: &Self) -> bool {
        (self + &(other * (F::ZERO - F::ONE)))
            .frobenius_norm()
            .epsilon_equals(&F::ZERO)
    }
}

#[rustfmt::skip::macros(vec, mat)]
#[cfg(test)]
mod tests {
    use super::Mat;
    use crate::core::num::{Field, Rational};

    macro_rules! mat {
		($($t:tt)*) => {
			mat_internal([$($t)*])
		};
	}

    fn mat_internal<F: Field, const R: usize, const C: usize>(arr: [[F; C]; R]) -> Mat<F> {
        Mat::new(R, C, arr.concat())
    }

    #[test]
    fn test_matrix_fmt_display() {
        let mat1 = mat![
			[0.0, 1.0, 2.0,  3.0],
			[4.0, 5.0, 6.0,  7.0],
			[8.0, 9.0, 10.0, 11.0]
		];

        println!("{}", mat1);
        assert_eq!(
            format!("{}", mat1),
            "\n[0.0, 1.0, 2.0, 3.0]\n[4.0, 5.0, 6.0, 7.0]\n[8.0, 9.0, 10.0, 11.0]"
        );
    }

    #[test]
    fn test_matrix_operations() {
        let mat1 = mat![
			[0.0, 1.0, 2.0,  3.0],
			[4.0, 5.0, 6.0,  7.0],
			[8.0, 9.0, 10.0, 11.0]
		];
        let mat2 = mat![
			[ 1.0, -1.0,  -2.0,  -3.0],
			[-4.0, -5.0,  -6.0,  -7.0],
			[-8.0, -9.0, -10.0, -11.0]
		];
        let mat3 = mat![
			[1.0, 0.0, 0.0, 0.0],
			[0.0, 0.0, 0.0, 0.0],
			[0.0, 0.0, 0.0, 0.0]
		];
        let mat4 = mat![
			[1.0, 0.0, 0.0],
			[0.0, 1.0, 0.0],
			[0.0, 0.0, 1.0]
		];
        let mat5 = mat![
			[4.0],
			[3.0],
			[2.0],
			[1.0],
		];
        let mat6 = mat![
			[10.0],
			[50.0],
			[90.0],
		];
        let mat7 = mat![
			[1.0],
			[5.0],
			[9.0],
		];

        assert_eq!(&mat1 + &mat2, mat3);
        assert_eq!(&mat4 * &mat1, mat1);
        assert_eq!(&mat4 * &mat2, mat2);

        assert_eq!(&mat1 * &mat5, mat6);
        assert_eq!(&mat7 * 10.0, mat6);
    }

    #[test]
    fn test_matrix_constructions() {
        let mat1 = Mat::new(
            3,
            4,
            vec![
				1.0, 0.0, 0.0, 0.0,
				0.0, 2.0, 0.0, 0.0,
				0.0, 0.0, 3.0, 0.0
			],
        );
        let mat2 = Mat::new(
            4,
            3,
            vec![
				1.0, 0.0, 0.0,
				0.0, 2.0, 0.0,
				0.0, 0.0, 3.0,
				0.0, 0.0, 0.0
			],
        );
        assert_eq!(mat1.transpose(), mat2);
        assert_eq!(mat1, mat2.transpose());

        let mat3 = Mat::new(
            4,
            4,
            vec![
				0.0, 1.0, 0.0, 0.0,
				1.0, 0.0, 0.0, 0.0,
				0.0, 0.0, 1.0, 0.0,
				0.0, 0.0, 0.0, 1.0
			],
        );
        assert_eq!(Mat::permutation(4, 0, 1), mat3);

        let mat4 = Mat::new(
            4,
            4,
            vec![
				1.0, 0.0, 0.0, 0.0,
				0.0, 1.0, 0.0, 0.0,
				2.0, 0.0, 1.0, 0.0,
				0.0, 0.0, 0.0, 1.0
			],
        );
        assert_eq!(Mat::replacement(4, 0, 2, 2.0), mat4);
    }

    #[test]
    fn test_matrix_operations_fractions() {
        let mat8 = Mat::new(
            2,
            2,
            vec![
				Rational::new_i64(1, 2), Rational::new_i64(1, 3),
				Rational::new_i64(1, 4), Rational::new_i64(1, 5)
			],
        );
        let mat9 = Mat::new(
            2,
            1,
            vec![
				Rational::new_i64(4, 1),
				Rational::new_i64(15, 1),
			],
        );
        let mat10 = Mat::new(
            2,
            1,
            vec![
				Rational::new_i64(7, 1),
				Rational::new_i64(4, 1),
			],
        );

        assert_eq!(&mat8 * &mat9, mat10);
    }
}
