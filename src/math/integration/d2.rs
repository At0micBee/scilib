use super::utils::*;

pub fn quad_fn<T>(
    function: impl Fn(f64, f64) -> f64,
    x_lower_bound: f64,
    x_upper_bound: f64,
    y_lower_bound: BoundaryD1<T>,
    y_upper_bound: BoundaryD1<T>,
    div: usize,
) where
    T: Fn(f64) -> f64,
{
}

pub fn testing() {
    quad_fn(
        |x, y| x * y,
        0.0,
        1.0,
        BoundaryD1::FixedValue(0.0),
        BoundaryD1::Function(|x| x),
        1000,
    );
}