//!
//! # Range creation tool
//! 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # Creating a range by increment
/// 
/// Creating a range with the start, stop, and increment. If the range is invalid, an empty Vec is returned.
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
