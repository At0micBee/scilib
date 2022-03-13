//!
//! # Range creation tool
//! 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # Creating a range by increment
/// 
/// Creating a range with the start, stop, and increment. If the range is invalid, an empty Vec is returned.
/// The points of this range are evenly spaced.
/// 
/// ```
/// # use scilib::range::by_increment;
/// let res = by_increment(2.0, -1.0, -0.5);
/// 
/// assert_eq!(res[0], 2.0);
/// assert_eq!(res.last().unwrap(), &-1.0);
/// assert_eq!(res.len(), 7);
/// 
/// ```
pub fn by_increment<T, U, V>(start: T, stop: U, step: V) -> Vec<f64>
where T: Into<f64>, U: Into<f64>, V: Into<f64> {

    let mut res: Vec<f64> = Vec::new();

    // Casting values
    let mut computed: f64 = start.into();
    let end: f64 = stop.into();
    let inc: f64 = step.into();

    // If the sign of this is pos, the increment doesn't make sense
    if ((computed - end) / inc).is_sign_positive() {
        return res;
    }

    'pushing: loop {
        if ((computed - end) / inc).is_sign_positive() {
            break 'pushing;
        }
        res.push(computed);
        computed += inc;
    }

    res
}

/// # Creating a range from a number of points
/// 
/// Creating a range with the start, stop, and number of points. If the range is invalid, an empty Vec is returned.
/// The points of this range are evenly spaced.
/// 
/// ```
/// # use scilib::range::linear;
/// let r = linear(0.0, -10.0, 15);
/// 
/// assert_eq!(r[0], 0.0);
/// assert_eq!(r.last().unwrap(), &-10.0);
/// assert_eq!(r.len(), 15);
/// ```
pub fn linear<T, U>(start: T, stop: U, n_points: usize) -> Vec<f64>
where T: Into<f64>, U: Into<f64> {
    
    let start_f: f64 = start.into();
    let stop_f: f64 = stop.into();
    let step: f64 = stop_f - start_f;
    let mut res: Vec<f64> = Vec::with_capacity(n_points);

    if n_points == 1 {
        res.push(start_f);
        return res;
    }

    for i in (0..n_points).map(|j| j as f64 / (n_points - 1) as f64) {
        res.push(start_f + i * step);
    }

    res
}

/// # Creating a logarithmic range from a number of points
/// 
/// Creating a range with the start, stop, and number of points. If the range is invalid, an empty Vec is returned.
/// The points of this range are space based on the log given log base.
/// 
/// ```
/// # use scilib::range::logarithmic;
/// let r = logarithmic(0.1, 10.0, 15, 10.0);
/// 
/// assert!((r[0] - 0.1).abs() < 1.0e-15);  // Start is correct
/// assert!((r[7] - 1.0).abs() < 1.0e-15);  // Distribution is correct
/// assert_eq!(r.last().unwrap(), &10.0);   // End is correct
/// assert_eq!(r.len(), 15);                // Number of points is correct
/// ```
pub fn logarithmic<T, U, V>(start: T, stop: U, n_points: usize, base: V) -> Vec<f64>
where T: Into<f64>, U: Into<f64>, V: Into<f64> {

    let start_f: f64 = start.into();
    let stop_f: f64 = stop.into();
    let base_f: f64 = base.into();

    // If the scale includes 0.0, it doesn't work
    if start_f == 0.0 || stop_f == 0.0 || start_f.signum() != stop_f.signum() {
        return vec![];
    }

    if n_points == 1 {
        return vec![start_f];
    }

    let sg: f64 = start_f.signum();
    let start_l: f64 = start_f.abs().log(base_f);
    let stop_l: f64 = stop_f.abs().log(base_f);

    let mut res: Vec<f64> = linear(start_l, stop_l, n_points);

    for elem in res.iter_mut() {
        *elem = sg * base_f.powf(*elem);
    }

    res
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # ∑ operation
/// 
/// Implementation of the mathematical summation.
/// * `type` - the type of the result
/// * `var` - name of the variable used in the loop (ex: i, n, j, etc)
/// * `range` - range used for the loop
/// * `form` - the formula summed
/// 
/// ```
/// # #[macro_use] extern crate scilib; 
/// # fn main() {
/// # use scilib::math::basic::*;
/// let res1 = summation!(u32, n, 0..=10, n);
/// assert_eq!(55, res1);
/// let res2 = summation!(f64, n, 0..=10, n as f64 * 0.5);
/// assert_eq!(27.5, res2);
/// # }
/// ```
#[macro_export]
macro_rules! summation {
    ($type:ty, $var:pat, $range:expr, $form:expr) => {{
        let mut result: $type = 0 as $type;
        for $var in $range {
            result += $form;
        }
        result
    }};
}

/// # Π operation
/// 
/// Implementation of the mathematical product.
/// * `type` - the type of the result
/// * `var` - name of the variable used in the loop (ex: i, n, j, etc)
/// * `range` - range used for the loop
/// * `form` - the formula multiplied
/// 
/// ```
/// # #[macro_use] extern crate scilib; 
/// # fn main() {
/// # use scilib::math::basic::*;
/// let res = product!(u32, i, 1..=10, i); // = 10!
/// assert_eq!(3628800, res);
/// # }
/// ```
#[macro_export]
macro_rules! product {
    ($type:ty, $var:pat, $range:expr, $form:expr) => {{
        let mut result: $type = 1 as $type;
        for $var in $range {
            result *= $form;
        }
        result
    }};
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
