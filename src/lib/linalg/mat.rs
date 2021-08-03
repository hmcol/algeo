use super::super::num::Field;
use super::util::get_box_iter;
use std::ops::{Add, Mul};
use std::fmt;

#[derive(Debug, PartialEq, Clone)]
pub struct Mat<F: Field> {
	entries: Vec<F>,
	rows: usize,
	cols: usize
}

impl<F: Field> Mat<F> {
	pub fn new(rows: usize, cols: usize, entries: Vec<F>) -> Mat<F> {
		Mat {
			entries,
			rows,
			cols
		}
	}

	pub fn identity(rows: usize) -> Mat<F> {
		Mat::new(rows, rows, 
			get_box_iter(rows, rows).map(|(r,c)|
				if r==c { F::ONE }
				else { F::ZERO }
			).collect()
		)
	}

	pub fn e(k: usize, n:usize) -> Mat<F> {
		Mat::new(n,1,
			(0..n).map(|r|
				if r==k { F::ONE }
				else { F::ZERO }
			).collect()
		)
	}

	pub fn zero(rows:usize, cols:usize) -> Mat<F> {
		Mat::new(rows, cols, 
			get_box_iter(rows, cols).map(|(_r,_c)|
				F::ZERO
			).collect()
		)
	}

	pub fn scale(rows: usize, i: usize, scale: F) -> Mat<F> {
		Mat::new(rows, rows, 
			get_box_iter(rows, rows).map(|(r,c)|
				if r==c && r==i { scale }
				else if r==c { F::ONE }
				else { F::ZERO }
			).collect()
		)
	}

	pub fn permutation(rows: usize, source: usize, target: usize) -> Mat<F> {
		Mat::new(rows, rows,
			get_box_iter(rows, rows).map(|(r,c)|
				if (r==source && c==target) || (r==target && c==source) { F::ONE }
				else if r==c && r != source && c != target { F::ONE }
				else { F::ZERO }
			).collect()
		)
	}

	pub fn replacement(rows: usize, source: usize, target: usize, scalar: F) -> Mat<F> {
		Mat::new(rows, rows,
			get_box_iter(rows, rows).map(|(r,c)|
				if c==source && r==target { scalar }
				else if r==c { F::ONE }
				else { F::ZERO }
			).collect()
		)
	}

	pub fn transpose(&self) -> Mat<F> {
		Mat::new(self.cols, self.rows,
			get_box_iter(self.cols, self.rows).map(|(c,r)| self.get(r,c).clone()).collect()
		)
	}

	pub fn get(&self, r: usize, c: usize) -> &F {
		&self.entries[r*self.cols+c]
	}

	pub fn entries(&self) -> &[F] {
		&self.entries
	}

	pub fn get_mut<'a>(&'a mut self, r: usize, c: usize) -> &'a mut F {
		&mut self.entries[r*self.cols+c]
	}

	pub fn rows(&self) -> usize { self.rows }
	pub fn cols(&self) -> usize { self.cols	}

	pub fn get_row_iter(&self, r: usize) -> impl Iterator<Item=&F> {
		self.entries[r*self.cols..(r+1)*self.cols].iter()
	}

	pub fn get_row_iter_mut(&mut self, r: usize) -> impl Iterator<Item=&mut F> {
		self.entries[r*self.cols..(r+1)*self.cols].iter_mut()
	}

	/// transposed to be column vector
	pub fn get_row(&self, r: usize) -> Mat<F> {
		Mat {
			entries: self.get_row_iter(r).map(|x| x.clone()).collect(),
			rows: self.cols,
			cols: 1
		}
	}

	pub fn get_col_iter(&self, c: usize) -> impl Iterator<Item=&F> {
		(0..self.rows).map(move |r| self.get(r,c))
	}

	pub fn get_col_iter_mut(& mut self, c: usize) -> impl Iterator<Item=&mut F> {
		// magic code that I got from discord.
		// the more simple (0..self.rows).map(move |r| entries.get_mut(r,c))
		// has weird lifetime issues. This works though, so sicko mode
		self.entries.chunks_exact_mut(self.rows).map(move |row| &mut row[c])
	}

	pub fn get_col(&self, c:usize) -> Mat<F> {
		Mat {
			entries: self.get_col_iter(c).map(|x| x.clone()).collect(),
			rows: self.rows,
			cols: 1
		}
	}
}

impl<F: Field> fmt::Display for Mat<F> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
		write!(f,"\n")?;
		for r in 0..self.rows {
			f.debug_list().entries(&self.entries[r*self.cols..r*self.cols+self.cols]).finish()?;
			if r!=self.rows-1 {
				write!(f,"\n")?;
			}
		}
		Ok(())
    }
}

impl<F: Field> Add for &Mat<F> {
	type Output = Mat<F>;

	fn add(self, other: Self) -> Mat<F>{
		if !(self.rows == other.rows && self.cols == other.cols) {
			panic!("tried to add matrix with incompatible dimensions.")
		}
		return Mat::new(self.rows, self.cols,
			self.entries.iter().zip(other.entries.iter()).map(|(&x,&y)| x+y).collect::<Vec<F>>()
		);
	}
}

// scalar multiplication
impl<F: Field> Mul<F> for &Mat<F> {
	type Output = Mat<F>;

	fn mul(self, scalar: F) -> Mat<F> {
		return Mat::new(self.rows, self.cols,
			self.entries.iter().map(|&x| x*scalar).collect()
		);
	}
}

impl<F: Field> Mul for &Mat<F> {
	type Output = Mat<F>;

	fn mul(self, other: Self) -> Mat<F>{
		if self.cols != other.rows {
			panic!("tried to multiply matrices with incompatible dimensions.")
		}
		Mat::new(self.rows, other.cols,
			get_box_iter(self.rows, other.cols)
				.map(|(r,c)|
					self.get_row_iter(r).zip(other.get_col_iter(c))
						.fold(F::ZERO,
							|accum: F, (&x,&y): (&F,&F)| accum + x*y) 
				).collect()
		)
	}
}








#[cfg(test)]
mod tests{
	use super::super::super::frac::Frac;
	use super::Mat;

	#[test]
	fn test_matrix_fmt_display(){
		let mat1 = Mat::new(3,4,
			vec![
				0.0, 1.0, 2.0,  3.0,
				4.0, 5.0, 6.0,  7.0,
				8.0, 9.0, 10.0, 11.0
			]
		);
	
		println!("{}",mat1);
		assert_eq!(
			format!("{}", mat1), 
			"\n[0.0, 1.0, 2.0, 3.0]\n[4.0, 5.0, 6.0, 7.0]\n[8.0, 9.0, 10.0, 11.0]"
		);
	}
	
	#[test]
	fn test_matrix_operations(){
		let mat1 = Mat::new(3,4,
			vec![
				0.0, 1.0, 2.0,  3.0,
				4.0, 5.0, 6.0,  7.0,
				8.0, 9.0, 10.0, 11.0
			]
		);
		let mat2 = Mat::new(3,4,
			vec![
				1.0, -1.0, -2.0,  -3.0,
				-4.0, -5.0, -6.0,  -7.0,
				-8.0, -9.0, -10.0, -11.0
			]
		);
		let mat3 = Mat::new(3,4,
			vec![
				1.0, 0.0, 0.0, 0.0,
				0.0, 0.0, 0.0, 0.0,
				0.0, 0.0, 0.0, 0.0
			]
		);
		let mat4 = Mat::new(3,3,
			vec![
				1.0, 0.0, 0.0,
				0.0, 1.0, 0.0,
				0.0, 0.0, 1.0
			]
		);
		let mat5 = Mat::new(4,1,
			vec![
				4.0,
				3.0,
				2.0,
				1.0,
			]
		);
		let mat6 = Mat::new(3,1,
			vec![
				10.0,
				50.0,
				90.0,
			]
		);
		let mat7 = Mat::new(3,1,
			vec![
				1.0,
				5.0,
				9.0,
			]
		);
	
		assert_eq!(&mat1+&mat2,mat3);
		assert_eq!(&mat4*&mat1, mat1);
		assert_eq!(&mat4*&mat2, mat2);
		
		assert_eq!(&mat1*&mat5, mat6);
		assert_eq!(&mat7*10.0, mat6);
	}
	
	#[test]
	fn test_matrix_constructions(){
		let mat1 = Mat::new(3,4,
			vec![
				1.0, 0.0, 0.0, 0.0,
				0.0, 2.0, 0.0, 0.0,
				0.0, 0.0, 3.0, 0.0
			]
		);
		let mat2 = Mat::new(4,3,
			vec![
				1.0, 0.0, 0.0,
				0.0, 2.0, 0.0,
				0.0, 0.0, 3.0,
				0.0, 0.0, 0.0
			]
		);
		assert_eq!(mat1.transpose(), mat2);
		assert_eq!(mat1, mat2.transpose());
	
	
		let mat3 = Mat::new(4,4,
			vec![
				0.0, 1.0, 0.0, 0.0,
				1.0, 0.0, 0.0, 0.0,
				0.0, 0.0, 1.0, 0.0,
				0.0, 0.0, 0.0, 1.0
			]
		);
		assert_eq!(Mat::permutation(4, 0, 1), mat3);
	
		let mat4 = Mat::new(4,4,
			vec![
				1.0, 0.0, 0.0, 0.0,
				0.0, 1.0, 0.0, 0.0,
				2.0, 0.0, 1.0, 0.0,
				0.0, 0.0, 0.0, 1.0
			]
		);
		assert_eq!(Mat::replacement(4, 0, 2, 2.0), mat4);
	}

	#[test]
	fn test_matrix_operations_fractions(){
		let mat8 = Mat::new(2,2,
			vec![
				Frac::new(1,2), Frac::new(1,3),
				Frac::new(1,4), Frac::new(1,5)
			]
		);
		let mat9 = Mat::new(2,1,
			vec![
				Frac::new(4,1),
				Frac::new(15,1),
			]
		);
		let mat10 = Mat::new(2,1,
			vec![
				Frac::new(7,1),
				Frac::new(4,1),
			]
		);

		assert_eq!(&mat8*&mat9, mat10);
	}
}

