use std::cmp::Ordering;

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