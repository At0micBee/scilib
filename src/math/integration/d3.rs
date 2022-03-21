pub fn quad_fn(
    function: impl Fn(f64, f64, f64) -> f64,
    x_lower_bound: f64,
    x_upper_bound: f64,
    y_lower_bound: impl Fn(f64) -> f64,
    y_upper_bound: impl Fn(f64) -> f64,
    z_lower_bound: impl Fn(f64, f64) -> f64,
    z_upper_bound: impl Fn(f64, f64) -> f64,
    div_x: usize,
    div_y: usize,
    div_z: usize,
) -> f64 {
    assert!(div_x >= 1 && div_y >= 1);
    let mut area = 0.0;
    let step_x = (x_upper_bound - x_lower_bound) / div_x as f64;
    for i in 0..div_x {
        let curr_x = x_lower_bound + i as f64 * step_x;
        let curr_y_lower_bound = y_lower_bound(curr_x);
        let step_y = (y_upper_bound(curr_x) - curr_y_lower_bound) / div_y as f64;
        let mut sub_area = 0.0;
        for j in 0..div_y {
            let curr_y = curr_y_lower_bound + j as f64 * step_y;
            let curr_z_lower_bound = z_lower_bound(curr_x, curr_y);
            let step_z = (z_upper_bound(curr_x, curr_y) - curr_z_lower_bound) / div_z as f64;
            let mut sub_sub_area = 0.0;
            for k in 0..div_z {
                sub_sub_area += function(curr_x, curr_y, curr_z_lower_bound + k as f64 * step_z);
            }
            sub_area += step_z * sub_sub_area;
        }
        area += step_y * sub_area;
    }
    step_x * area
}
