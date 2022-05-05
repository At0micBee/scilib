////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # (d3) Rectangle rule
/// Integrate three-dimensional function through the rectangle rule:
/// * `function` - Closure or function pointer matching `f(x, y, z) = w`
/// * `x_lower_bound` - lower bound of the first integral (Fixed value)
/// * `x_upper_bound` - upper bound of the first integral (Fixed value)
/// * `y_lower_bound` - lower bound of the second integral (Closure or function pointer `f(x) = y`)
/// * `y_upper_bound` - upper bound of the second integral (Closure or function pointer `f(x) = y`)
/// * `z_lower_bound` - lower bound of the third integral (Closure or function pointer `f(x, y) = z`)
/// * `z_upper_bound` - upper bound of the third integral (Closure or function pointer `f(x, y) = z`)
/// * `div_x` - Number of chunk evaluated for the first variable: _big_ = high precision & low performances, _small_ = high performances & low precision
/// * `div_y` - Same as `div_x` for the second variable
/// * `div_z` - Same as `div_x` for the third variable
///
/// ```
/// # use scilib::math::integration::d3::*;
/// let res = rect(|x, y, z| x * y * z, -3.0, 5.0, |x| x, |x| 2.0 * x, |_x, y| y, |_x, y| 2.0 * y, 10000, 1000, 100);
/// assert!((res - 13965.0).abs() < 70.0);
/// ```
pub fn rect(
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

fn trapz_sub_sub_z(
    function: impl Fn(f64, f64, f64) -> f64 + Copy,
    curr_x: f64,
    curr_y: f64,
    z_lower_bound: impl Fn(f64, f64) -> f64 + Copy,
    z_upper_bound: impl Fn(f64, f64) -> f64 + Copy,
    div_z: usize,
) -> f64 {
    let curr_z_lower_bound = z_lower_bound(curr_x, curr_y);
    let step = (z_upper_bound(curr_x, curr_y) - curr_z_lower_bound) / div_z as f64;
    let mut area = 0.0;
    let mut save_previous = function(curr_x, curr_y, curr_z_lower_bound);
    for k in 0..div_z {
        let next_value = function(curr_x, curr_y, curr_z_lower_bound + (k + 1) as f64 * step);
        area += save_previous + next_value;
        save_previous = next_value;
    }
    step * area / 2.0
}

fn trapz_sub_y(
    function: impl Fn(f64, f64, f64) -> f64 + Copy,
    curr_x: f64,
    y_lower_bound: impl Fn(f64) -> f64 + Copy,
    y_upper_bound: impl Fn(f64) -> f64 + Copy,
    z_lower_bound: impl Fn(f64, f64) -> f64 + Copy,
    z_upper_bound: impl Fn(f64, f64) -> f64 + Copy,
    div_y: usize,
    div_z: usize,
) -> f64 {
    let curr_y_lower_bound = y_lower_bound(curr_x);
    let step = (y_upper_bound(curr_x) - curr_y_lower_bound) / div_y as f64;
    let mut area = 0.0;
    let mut save_previous = trapz_sub_sub_z(
        function,
        curr_x,
        curr_y_lower_bound,
        z_lower_bound,
        z_upper_bound,
        div_z,
    );
    for j in 0..div_y {
        let next_value = trapz_sub_sub_z(
            function,
            curr_x,
            curr_y_lower_bound + (j + 1) as f64 * step,
            z_lower_bound,
            z_upper_bound,
            div_z,
        );
        area += save_previous + next_value;
        save_previous = next_value;
    }
    step * area / 2.0
}

/// # (d3) Trapezoidal rule
/// Integrate three-dimensional function through the trapezoidal rule:
/// * `function` - Closure or function pointer matching `f(x, y, z) = w`
/// * `x_lower_bound` - lower bound of the first integral (Fixed value)
/// * `x_upper_bound` - upper bound of the first integral (Fixed value)
/// * `y_lower_bound` - lower bound of the second integral (Closure or function pointer `f(x) = y`)
/// * `y_upper_bound` - upper bound of the second integral (Closure or function pointer `f(x) = y`)
/// * `z_lower_bound` - lower bound of the third integral (Closure or function pointer `f(x, y) = z`)
/// * `z_upper_bound` - upper bound of the third integral (Closure or function pointer `f(x, y) = z`)
/// * `div_x` - Number of chunk evaluated for the first variable: _big_ = high precision & low performances, _small_ = high performances & low precision
/// * `div_y` - Same as `div_x` for the second variable
/// * `div_z` - Same as `div_x` for the third variable
///
/// ```
/// # use scilib::math::integration::d3::*;
/// let res = trapz(|x, y, z| x * y * z, -3.0, 5.0, |x| x, |x| 2.0 * x, |_x, y| y, |_x, y| 2.0 * y, 10000, 1000, 100);
/// assert!((res - 13965.0).abs() < 10e-2);
/// ```
pub fn trapz(
    function: impl Fn(f64, f64, f64) -> f64 + Copy,
    x_lower_bound: f64,
    x_upper_bound: f64,
    y_lower_bound: impl Fn(f64) -> f64 + Copy,
    y_upper_bound: impl Fn(f64) -> f64 + Copy,
    z_lower_bound: impl Fn(f64, f64) -> f64 + Copy,
    z_upper_bound: impl Fn(f64, f64) -> f64 + Copy,
    div_x: usize,
    div_y: usize,
    div_z: usize,
) -> f64 {
    let step = (x_upper_bound - x_lower_bound) / div_x as f64;
    let mut area = 0.0;
    let mut save_previous = trapz_sub_y(
        function,
        x_lower_bound,
        y_lower_bound,
        y_upper_bound,
        z_lower_bound,
        z_upper_bound,
        div_y,
        div_z,
    );
    for i in 0..div_x {
        let next_value = trapz_sub_y(
            function,
            x_lower_bound + (i + 1) as f64 * step,
            y_lower_bound,
            y_upper_bound,
            z_lower_bound,
            z_upper_bound,
            div_y,
            div_z,
        );
        area += save_previous + next_value;
        save_previous = next_value;
    }
    step * area / 2.0
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
