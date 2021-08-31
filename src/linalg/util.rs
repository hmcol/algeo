use std::cmp::Ordering;

use crate::core::num::Field;

use super::mat::Mat;

// (0,0), (0,1), ..., (0,cols)
// (1,0), (1,1), ..., (1,cols)
//  ...                 ...
// (row,0) ...
pub fn get_box_iter(rows: usize, cols: usize) -> impl Iterator<Item=(usize, usize)>{
	(0..rows).map(move |r| ((0..cols).map(move |c| (r,c) ) ) ).flatten()
}

pub fn get_max_index<T, I, C>(slice: I, le: C) -> usize where
	I: Iterator<Item=T>,
	C: Fn(&T, &T)->Option<Ordering> {

	slice
		.enumerate()
        .max_by(|(_n, x), (_m, y)| le(x, y).unwrap())
        .unwrap()
		.0
}

pub fn num_to_seq(n: usize, base: usize) -> impl Iterator<Item=(usize)> {
	BaseSeqIterator::new(n, base)
}

struct BaseSeqIterator {
	n: usize,
	base: usize
}

impl BaseSeqIterator {
	pub fn new(n: usize, base: usize) -> BaseSeqIterator {
		BaseSeqIterator {
			n: n*base,
			base
		}
	}
}

impl Iterator for BaseSeqIterator {
	type Item = usize;

	fn next(&mut self) -> Option<Self::Item> {
		self.n /= self.base;

		if self.n == 0 {
			None
		} else {
			Some(self.n % self.base)
		}
	}
}

pub fn mat_iterator<'a, F: Field>(n: usize, m: usize, values: &'a [F])-> impl 'a + Iterator<Item=Mat<F>>{
	MatIterator {
		values,
		n,
		m,
		num: 0
	}
}

pub struct MatIterator<'a, F: Field> {
	values: &'a [F],
	n: usize,
	m: usize,
	num: usize,
}

impl<'a, F: Field> Iterator for MatIterator<'a, F> {
    type Item = Mat<F>;

    fn next(&mut self) -> Option<Mat<F>> {
        if self.num >= self.values.len().pow((self.n * self.m) as u32) {
			None
		} else {
			let mut entries: Vec<F> = num_to_seq(self.num, self.values.len())
				.map(|i| self.values[i]).collect();
			
			for _i in (entries.len())..(self.n * self.m) {
				entries.push(self.values[0]);
			}

			self.num += 1;

			Some(Mat::new(self.n, self.m, entries))
		}
    }
}