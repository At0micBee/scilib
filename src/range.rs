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
pub fn by_increment(start: f64, stop: f64, inc: f64) -> Vec<f64> {

    let mut res: Vec<f64> = Vec::new();

    // If the sign of this is pos, the increment doesn't make sense
    if ((start - stop) / inc).is_sign_positive() {
        return res;
    }

    let mut computed: f64 = start;

    'pushing: loop {
        if ((computed - stop) / inc).is_sign_positive() {
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
pub fn linear(start: f64, stop: f64, n_points: usize) -> Vec<f64> {
    
    let step: f64 = stop - start;
    let mut res: Vec<f64> = Vec::with_capacity(n_points);

    if n_points == 1 {
        res.push(start);
        return res;
    }

    for i in (0..n_points).map(|j| j as f64 / (n_points - 1) as f64) {
        res.push(start + i * step);
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
pub fn logarithmic(start: f64, stop: f64, n_points: usize, base: f64) -> Vec<f64> {

    // If the scale includes 0.0, it doesn't work
    if start == 0.0 || stop == 0.0 || start.signum() != stop.signum() {
        return vec![];
    }

    if n_points == 1 {
        return vec![start];
    }

    let sg: f64 = start.signum();
    let start_l: f64 = start.abs().log(base);
    let stop_l: f64 = stop.abs().log(base);

    let mut res: Vec<f64> = linear(start_l, stop_l, n_points);

    for elem in res.iter_mut() {
        *elem = sg * base.powf(*elem);
    }

    res
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
