#![allow(unused)]

use std::{ops::{Index, IndexMut}, fmt};

/// internal type for storing 2d data
/// 
/// different from like mathematical matrix operations
#[derive(Debug, PartialEq, Clone)]
pub struct Mat<T> {
    entries: Vec<T>,
    rows: usize,
    cols: usize,
}


impl<T> Mat<T> {
    pub fn new(rows: usize, cols: usize, entries: Vec<T>) -> Mat<T> {
        Mat {
            entries,
            rows,
            cols,
        }
    }

    pub fn swap(&mut self, (r1, c1): (usize, usize), (r2, c2): (usize, usize)) {
        self.entries.swap(r1 * self.cols + c1, r2 * self.cols + c2);
    }

    // Mat<T> constructors ----------------------------------------------------

    /// identity matrix with of dimension `n * n`
    // pub fn identity(n: usize) -> Mat<F> {
    //     Mat::new(
    //         n,
    //         n,
    //         get_box_iter(n, n)
    //             .map(|(r, c)| if r == c { F::ONE } else { F::ZERO })
    //             .collect(),
    //     )
    // }

    

    /// `k`th standard basis vector of dimension `n`
    // pub fn e(k: usize, n: usize) -> Mat<T> {
    //     Mat::new(
    //         n,
    //         1,
    //         (0..n)
    //             .map(|r| if r == k { F::ONE } else { F::ZERO })
    //             .collect(),
    //     )
    // }

    /// zero matrix of dimension `rows * cols`
    // pub fn zero(rows: usize, cols: usize) -> Mat<F> {
    //     Mat::new(
    //         rows,
    //         cols,
    //         get_box_iter(rows, cols).map(|(_r, _c)| F::ZERO).collect(),
    //     )
    // }

    /// elementary matrix of dimension `rows*rows` that scales the `i`th row
    // pub fn scale_mat(rows: usize, i: usize, scale: F) -> Mat<F> {
    //     Mat::new(
    //         rows,
    //         rows,
    //         get_box_iter(rows, rows)
    //             .map(|(r, c)| {
    //                 if r == c && r == i {
    //                     scale
    //                 } else if r == c {
    //                     F::ONE
    //                 } else {
    //                     F::ZERO
    //                 }
    //             })
    //             .collect(),
    //     )
    // }

    /// elementary matrix of dimension `rows*rows` that permutes the `source`
    /// and `target` rows.
    // pub fn permutation(rows: usize, source: usize, target: usize) -> Mat<F> {
    //     Mat::new(
    //         rows,
    //         rows,
    //         get_box_iter(rows, rows)
    //             .map(|(r, c)| {
    //                 if (r, c) == (source, target)
    //                     || (c, r) == (source, target)
    //                     || (r == c && r != source && c != target)
    //                 {
    //                     F::ONE
    //                 } else {
    //                     F::ZERO
    //                 }
    //             })
    //             .collect(),
    //     )
    // }

    /// input - permutation : Vec<F> of length n with entries `0..(n-1)` in
    /// some order. The vec can be thought of as a bijection of
    /// {0,1,...,n-1} -> {0,1,...,n-1}
    /// and `permutation_from_vec` returns a matrix that does this permutation
    /// to the rows.
    // pub fn permutation_from_vec(permutation: Vec<usize>) -> Mat<F> {
    //     let n = permutation.len();
    //     Mat::new(
    //         n,
    //         n,
    //         get_box_iter(n, n)
    //             .map(|(r, c)| if permutation[r] == c { F::ONE } else { F::ZERO })
    //             .collect(),
    //     )
    // }

    /// elementary matrix of dimension `rows*rows` adds the `source` row
    /// scaled by `scalar` to the `target` row.
    // pub fn replacement(rows: usize, source: usize, target: usize, scalar: T) -> Mat<T> {
    //     Mat::new(
    //         rows,
    //         rows,
    //         get_box_iter(rows, rows)
    //             .map(|(r, c)| {
    //                 if c == source && r == target {
    //                     scalar
    //                 } else if r == c {
    //                     F::ONE
    //                 } else {
    //                     F::ZERO
    //                 }
    //             })
    //             .collect(),
    //     )
    // }

    // Mat<F> mutators ---------------------------------------------------------

    // pub fn scale_row(&mut self, i: usize, scalar: F) {
    //     for c in 0..self.cols {
    //         self[(i, c)] = scalar * self[(i, c)];
    //     }
    // }

    

    // pub fn replace_row(&mut self, source: usize, target: usize, scalar: F) {
    //     for c in 0..self.cols {
    //         self[(target, c)] = self[(target, c)] + self[(source, c)] * scalar;
    //     }
    // }

    // TODO Refactor into MulAssign
    // pub fn scale(&mut self, scalar: F) {
    //     self.entries.iter_mut().for_each(|x| *x = *x * scalar);
    // }

    // Mat<F> getters (possibly mutable) ---------------------------------------

    pub fn get_unchecked(&self, r: usize, c: usize) -> &T {
        &self.entries[r * self.cols + c]
    }

    /// should replace with custom iterator
    pub fn entries(&self) -> &Vec<T> {
        &self.entries
    }

    pub fn get_mut_unchecked(&mut self, r: usize, c: usize) -> &mut T {
        &mut self.entries[r * self.cols + c]
    }

    pub fn rows(&self) -> usize {
        self.rows
    }
    pub fn cols(&self) -> usize {
        self.cols
    }

    pub fn size(&self) -> (usize, usize) {
        (self.rows, self.cols)
    }

    /// should replace with custom iterator
    pub fn get_row_iter(&self, r: usize) -> impl Iterator<Item = &T> {
        self.entries[r * self.cols..(r + 1) * self.cols].iter()
    }

    /// should replace with custom iterator
    pub fn get_row_iter_mut(&mut self, r: usize) -> impl Iterator<Item = &mut T> {
        self.entries[r * self.cols..(r + 1) * self.cols].iter_mut()
    }

    

    /// should replace with custom iterator
    pub fn get_col_iter(&self, c: usize) -> impl Iterator<Item = &T> {
        (0..self.rows).map(move |r| self.get_unchecked(r, c))
    }

    /// should replace with custom iterator
    pub fn get_col_iter_mut(&mut self, c: usize) -> impl Iterator<Item = &mut T> {
        // magic code that I got from discord.
        // the more simple (0..self.rows).map(move |r| entries.get_mut(r,c))
        // has weird lifetime issues. This works though, so sicko mode
        self.entries
            .chunks_exact_mut(self.rows)
            .map(move |row| &mut row[c])
    }

    

    // pub fn frobenius_norm(&self) -> F {
    //     self.entries.iter().fold(F::ZERO, |accum, &x| accum + x * x)
    // }

    // Operations --------------------------------------------------------------

    // pub fn dot(&self, other: &Mat<T>) -> T {
    //     self.entries
    //         .iter()
    //         .zip(other.entries())
    //         .fold(F::ZERO, |accum, (&x, &y)| accum + x * y)
    // }

    // Util --------------------------------------------------------------------

    // pub fn index_iter(&self) -> impl Iterator<Item = (usize, usize)> {
    //     (0..self.rows()).cartesian_product(0..self.cols())
    // }

    // Debug/Testing methods ---------------------------------------------------

    // / TODO: decide whether `is_upper_triangular` should require F to implement
    // / `EpsilonEquality`. In practice, I'm only using `is_upper_triangular` to
    // / verify the correctness of the LU decomposition, and in that case the
    // / replacement operations are exact (so that I wouldn't have to loosen
    // / equality to epsilon equality).
    // pub fn is_upper_triangular(&self) -> bool {
    //     for r in 0..self.rows {
    //         for c in 0..(std::cmp::min(r, self.cols)) {
    //             if self[(r, c)] != F::ZERO {
    //                 return false;
    //             }
    //         }
    //     }

    //     true
    // }

    // pub fn is_lower_triangular(&self) -> bool {
    //     for r in 0..self.rows {
    //         for c in (std::cmp::min(r, self.cols) + 1)..self.cols {
    //             if self[(r, c)] != F::ZERO {
    //                 return false;
    //             }
    //         }
    //     }

    //     true
    // }
}


impl<T: Copy> Mat<T> {
    pub fn from_col_vectors(vs: &[Mat<T>]) -> Mat<T> {
        Mat::new(
            vs[0].rows(),
            vs.len(),
            get_box_iter(vs[0].rows(), vs.len())
                .map(|(r, c)| vs[c][(r, 0)])
                .collect(),
        )
    }

    pub fn from_row_vectors(vs: &[Mat<T>]) -> Mat<T> {
        Mat::new(
            vs.len(),
            vs[0].rows(),
            get_box_iter(vs.len(), vs[0].rows())
                .map(|(r, c)| vs[r][(c, 0)])
                .collect(),
        )
    }


    pub fn permute_rows(&mut self, source: usize, target: usize) {
        for c in 0..self.cols {
            let temp = self[(source, c)];
            self[(source, c)] = self[(target, c)];
            self[(target, c)] = temp;
        }
    }
}

impl<T: Clone> Mat<T> {
    

    /// returns transposition of the matrix.
    pub fn transpose(&self) -> Mat<T> {
        Mat::new(
            self.cols,
            self.rows,
            get_box_iter(self.cols, self.rows)
                .map(|(c, r)| self.get_unchecked(r, c))
                .cloned()
                .collect(),
        )
    }

    /// transposed to be column vector
    pub fn get_row(&self, r: usize) -> Mat<T> {
        Mat::from(self.get_row_iter(r).cloned())
    }




    pub fn get_col(&self, c: usize) -> Mat<T> {
        Mat::from(self.get_col_iter(c).cloned())
    }
}

impl<T, I: Iterator<Item = T>> From<I> for Mat<T> {
    fn from(iter: I) -> Self {
        let entries: Vec<T> = iter.collect();
        let rows = entries.len();
        Mat {
            entries,
            rows,
            cols: 1,
        }
    }
}

impl<T> Index<(usize, usize)> for Mat<T> {
    type Output = T;

    fn index(&self, (r, c): (usize, usize)) -> &T {
        self.get_unchecked(r, c)
    }
}

impl<T> IndexMut<(usize, usize)> for Mat<T> {
    fn index_mut(&mut self, (r, c): (usize, usize)) -> &mut T {
        self.get_mut_unchecked(r, c)
    }
}

impl<T: fmt::Debug> fmt::Display for Mat<T> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f)?;
        for r in 0..self.rows {
            f.debug_list()
                .entries(self.get_row_iter(r))
                .finish()?;
            if r != self.rows - 1 {
                writeln!(f)?;
            }
        }
        Ok(())
    }
}


pub fn get_box_iter(rows: usize, cols: usize) -> impl Iterator<Item = (usize, usize)> {
    (0..rows).flat_map(move |r| ((0..cols).map(move |c| (r, c))))
}