pub fn quad_fn(
    function: impl Fn(f64, f64) -> f64,
    x_lower_bound: f64,
    x_upper_bound: f64,
    y_lower_bound: impl Fn(f64) -> f64,
    y_upper_bound: impl Fn(f64) -> f64,
    div_x: usize,
    div_y: usize,
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
            sub_area += function(curr_x, curr_y_lower_bound + j as f64 * step_y);
        }
        area += step_y * sub_area;
    }
    step_x * area
}

fn trapz_fn_sub_y(
    function: impl Fn(f64, f64) -> f64,
    curr_x: f64,
    y_lower_bound: impl Fn(f64) -> f64 + Copy,
    y_upper_bound: impl Fn(f64) -> f64 + Copy,
    div_y: usize,
) -> f64 {
    let curr_y_lower_bound = y_lower_bound(curr_x);
    let step = (y_upper_bound(curr_x) - curr_y_lower_bound) / div_y as f64;
    let mut save_previous = function(curr_x, curr_y_lower_bound);
    let mut area = 0.0;
    for j in 0..div_y {
        let next_value = function(curr_x, curr_y_lower_bound + (j + 1) as f64 * step);
        area += save_previous + next_value;
        save_previous = next_value;
    }
    step * area / 2.0
}

pub fn trapz_fn(
    function: impl Fn(f64, f64) -> f64 + Copy,
    x_lower_bound: f64,
    x_upper_bound: f64,
    y_lower_bound: impl Fn(f64) -> f64 + Copy,
    y_upper_bound: impl Fn(f64) -> f64 + Copy,
    div_x: usize,
    div_y: usize,
) -> f64 {
    assert!(div_x >= 1 && div_y >= 1);

    let step = (x_upper_bound - x_lower_bound) / div_x as f64;
    let mut save_previous =
        trapz_fn_sub_y(function, x_lower_bound, y_lower_bound, y_upper_bound, div_y);
    let mut area = 0.0;

    for i in 0..div_x {
        let next_value = trapz_fn_sub_y(
            function,
            x_lower_bound + (i + 1) as f64 * step,
            y_lower_bound,
            y_upper_bound,
            div_y,
        );
        area += next_value + save_previous;
        save_previous = next_value;
    }
    step * area / 2.0
}
