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
    y_lower_bound: f64,
    y_upper_bound: f64,
    div_y: usize,
) -> f64 {
    let step_y = (y_upper_bound - y_lower_bound) / div_y as f64;
    let mut save_y = function(curr_x, y_lower_bound);
    let mut sub_area = 0.0;
    for j in 0..div_y {
        let next_value_y = function(curr_x, y_lower_bound + (j + 1) as f64 * step_y);
        sub_area += save_y + next_value_y;
        save_y = next_value_y;
    }
    sub_area
}

pub fn trapz_fn(
    function: impl Fn(f64, f64) -> f64 + Copy,
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
    let mut save_x = trapz_fn_sub_y(
        function,
        x_lower_bound,
        y_lower_bound(x_lower_bound),
        y_upper_bound(x_lower_bound),
        div_y,
    );

    for i in 0..div_x {
        let curr_x = x_lower_bound + i as f64 * step_x;
        let next_x = curr_x + step_x;
        let curr_y_lower_bound = y_lower_bound(curr_x);
        let step_y = (y_upper_bound(curr_x) - curr_y_lower_bound) / div_y as f64;
        let next_value_x = trapz_fn_sub_y(
            function,
            next_x,
            y_lower_bound(next_x),
            y_upper_bound(next_x),
            div_y,
        );
        area += step_y * (save_x + next_value_x) / 2.0;
        save_x = next_value_x;
    }
    step_x * area / 2.0
}
