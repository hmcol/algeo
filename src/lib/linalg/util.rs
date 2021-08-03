
// (0,0), (0,1), ..., (0,cols)
// (1,0), (1,1), ..., (1,cols)
//  ...                 ...
// (row,0) ...
pub fn get_box_iter(rows: usize, cols: usize) -> impl Iterator<Item=(usize, usize)>{
	(0..rows).map(move |r| ((0..cols).map(move |c| (r,c) ) ) ).flatten()
}

pub fn get_max_index(slice: impl Iterator<Item=f64>) -> usize {
	// slice.enumerate().fold((0,f32::MIN),
	// 	|(max_index, max_val), (index,val)|
	// 		if val > max_val {
	// 			(index, val)
	// 		} else {
	// 			(max_index, max_val)
	// 		}
	// ).0;

	slice
		.enumerate()
        .max_by(|(x,_n), (y,_m)| x.partial_cmp(&y).unwrap())
        .unwrap()
		.0
}