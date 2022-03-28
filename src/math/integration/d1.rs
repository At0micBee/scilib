
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

fn chunk_rect(a: f64, offset: f64) -> f64 {
    a * offset
}

fn chunk_trapz(a: f64, b: f64, offset: f64) -> f64 {
    (a + b) / 2.0 * offset
}

fn chunk_simp(
    p0: (f64 /*a0*/, f64 /*b0*/),
    p1: f64, /*b1*/
    p2: (f64 /*a2*/, f64 /*b2*/),
) -> f64 {
    (p2.0 - p0.0) / 6.0 * (p0.1 + 4.0 * p1 + p2.1)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # (d1) Rectangle rule - Function
/// Integrate one-dimensional function through the rectangle rule:
/// * `function` - Closure or function pointer matching `f(x) = y`
/// * `lower_bound` - First limit of the integral (Fixed value)
/// * `upper_bound` - Second limit of the integral (Fixed value)
/// * `div` - Number of chunk evaluated: _big_ = high precision & low performances, _small_ = high performances & low precision
/// 
/// ```
/// # use scilib::math::integration::d1::*;
/// # use std::f64::consts::PI;
/// let res = rect_fn(|x| x.cos(), PI * 2.0, PI * 2.0 + PI / 2.0, 1000);
/// assert!((res - 1.0).abs() < 1.0e-3);
/// ```
pub fn rect_fn(
    function: impl Fn(f64) -> f64,
    lower_bound: f64,
    upper_bound: f64,
    div: usize,
) -> f64 {
    assert!(div >= 1);
    let mut area = 0.0;
    let step = (upper_bound - lower_bound) / div as f64;
    for i in 0..div {
        area += chunk_rect(function(step * (i as f64) + lower_bound), step);
    }
    area
}

/// # (d1) Trapezoidal rule - Function
/// Integrate one-dimensional function through the trapezoidal rule:
/// * `function` - Closure or function pointer matching `f(x) = y`
/// * `lower_bound` - First limit of the integral (Fixed value)
/// * `upper_bound` - Second limit of the integral (Fixed value)
/// * `div` - Number of chunk evaluated: _big_ = high precision & low performances, _small_ = high performances & low precision
/// 
/// ```
/// # use scilib::math::integration::d1::*;
/// let res = trapz_fn(|x| x.exp(), -2.0, 10.0, 10_000);
/// assert!((res - 22026.3304).abs() < 1.0e-2);
/// ```
pub fn trapz_fn(
    function: impl Fn(f64) -> f64,
    lower_bound: f64,
    upper_bound: f64,
    div: usize,
) -> f64 {
    assert!(div >= 1);
    let mut area = 0.0;
    let step = (upper_bound - lower_bound) / div as f64;
    let mut a = function(lower_bound);
    for i in 0..div {
        let b = function(step * (i as f64 + 1.0) + lower_bound);
        area += chunk_trapz(a, b, step);
        a = b;
    }
    area
}

/// # (d1) Simpson's rule - Function
/// Integrate one-dimensional function through the simpson rule:
/// * `function` - Closure or function pointer matching `f(x) = y`
/// * `lower_bound` - First limit of the integral (Fixed value)
/// * `upper_bound` - Second limit of the integral (Fixed value)
/// * `div` - Number of chunk evaluated: _big_ = high precision & low performances, _small_ = high performances & low precision
///
/// ```
/// # use scilib::math::integration::d1::*;
/// let res = simp_fn(|x| x.exp(), -2.0, 10.0, 1000);
/// assert!((res - 22026.330459).abs() < 1.0e-4);
/// ```
pub fn simp_fn(
    function: impl Fn(f64) -> f64,
    lower_bound: f64,
    upper_bound: f64,
    div: usize,
) -> f64 {
    assert!(div >= 2 && div % 2 == 0);
    let mut area = 0.0;
    let step = (upper_bound - lower_bound) / div as f64;
    let mut p0 = (lower_bound, function(lower_bound));
    for i in (0..div).step_by(2) {
        let curr_pos = lower_bound + (i + 2) as f64 * step;
        let p2 = (curr_pos, function(curr_pos));
        area += chunk_simp(p0, function((p0.0 + p2.0) / 2.0), p2);
        p0 = p2;
    }
    area
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # (d1) Rectangle rule - Data by tuples
/// Integrate one-dimensional function represented by data enclosed in tuples through the rectangle rule:
/// * `data` - array of tuples corresponding to `x` and `y` rescpectively
/// 
/// ```
/// # use scilib::math::integration::d1::*;
/// let count = 1_000;
/// let mut data = vec![(0.0, 0.0); count];
/// let lower_bound = -10.0;
/// let upper_bound = 10.0;
/// let step = (upper_bound - lower_bound) / (count - 1) as f64;
/// for i in 0..count {
///     let x = lower_bound + i as f64 * step;
///     data[i] = (x, x.sin());
/// }
/// let res = rect_dt_tup(&data);
/// assert!(res.abs() < 1.0e-1);
/// ```
pub fn rect_dt_tup(data: &[(f64, f64)]) -> f64 {
    let mut area = 0.0;
    for i in 0..data.len() - 1 {
        area += chunk_rect(data[i].1, data[i + 1].0 - data[i].0);
    }
    area
}

/// # (d1) Trapezoidal rule - Data by tuples
/// Integrate one-dimensional function represented by data enclosed in tuples through the trapezoidal rule:
/// * `data` - array of tuples corresponding to `x` and `y` rescpectively
/// 
/// ```
/// # use scilib::math::integration::d1::*;
/// let count = 100;
/// let mut data = vec![(0.0, 0.0); count];
/// let lower_bound = -10.0;
/// let upper_bound = 10.0;
/// let step = (upper_bound - lower_bound) / (count - 1) as f64;
/// for i in 0..count {
///     let x = lower_bound + i as f64 * step;
///     data[i] = (x, x.sin());
/// }
/// let res = trapz_dt_tup(&data);
/// assert!(res.abs() < 1.0e-3);
/// ```
pub fn trapz_dt_tup(data: &[(f64, f64)]) -> f64 {
    let mut area = 0.0;
    for i in 0..data.len() - 1 {
        area += chunk_trapz(data[i].1, data[i + 1].1, data[i + 1].0 - data[i].0);
    }
    area
}

fn chunk_simp_dt(
    p0: (f64 /*a0*/, f64 /*b0*/),
    p1: (f64 /*a1*/, f64 /*b1*/),
    p2: (f64 /*a2*/, f64 /*b2*/),
) -> f64 {
    let h0 = p1.0 - p0.0;
    let h1 = p2.0 - p1.0;
    (h0 + h1) / 6.0 * ((2.0 - h1 / h0) * p0.1 + (h0 + h1).powi(2) / (h0 * h1) * p1.1 + (2.0 - h0 / h1) * p2.1)
}

/// # (d1) Simpson's rule - Data by tuples
/// Integrate one-dimensional function represented by data enclosed in tuples through the trapezoidal rule:
/// * `data` - array of tuples corresponding to `x` and `y` rescpectively
/// 
/// **CAUTION**: `data.len()` must be an uneven number greater than 3
/// 
/// ```
/// # use scilib::math::integration::d1::*;
/// let count = 11;
/// let mut data = vec![(0.0, 0.0); count];
/// let lower_bound = -10.0;
/// let upper_bound = 10.0;
/// let step = (upper_bound - lower_bound) / (count - 1) as f64;
/// for i in 0..count {
///     let x = lower_bound + i as f64 * step;
///     data[i] = (x, x.sin());
/// }
/// let res = simp_dt_tup(&data);
/// assert!(res.abs() < 1.0e-4);
/// ```
pub fn simp_dt_tup(data: &[(f64, f64)]) -> f64 {
    assert!(data.len() >= 3 && data.len() % 2 == 1);
    let mut area = 0.0;
    for i in (0..data.len() - 2).step_by(2) {
        area += chunk_simp_dt(data[i], data[i + 1], data[i + 2]);
    }
    area
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # (d1) Rectangle rule - Data by slices
/// Integrate one-dimensional function represented by data enclosed in slices through the rectangle rule:
/// * `x` - array of x-coordinates
/// * `y` - array of y-coordinates
/// 
/// ```
/// # use scilib::math::integration::d1::*;
/// let count = 1_000;
/// let mut x = vec![0.0; count];
/// let mut y = vec![0.0; count];
/// let lower_bound = -10.0;
/// let upper_bound = 10.0;
/// let step = (upper_bound - lower_bound) / (count - 1) as f64;
/// for i in 0..count {
///     let pos = lower_bound + i as f64 * step;
///     x[i] = pos;
///     y[i] = pos.sin();
/// }
/// let res = rect_dt_sli(&x, &y);
/// assert!(res.abs() < 1.0e-1);
/// ```
pub fn rect_dt_sli(x: &[f64], y: &[f64]) -> f64 {
    assert!(x.len() == y.len());
    let mut area = 0.0;
    for i in 0..x.len() - 1 {
        area += chunk_rect(y[i], x[i + 1] - x[i]);
    }
    area
}

/// # (d1) Trapezoidal rule - Data by slices
/// Integrate one-dimensional function represented by data enclosed in slices through the trapezoidal rule:
/// * `x` - array of x-coordinates
/// * `y` - array of y-coordinates
/// 
/// ```
/// # use scilib::math::integration::d1::*;
/// let count = 100;
/// let mut x = vec![0.0; count];
/// let mut y = vec![0.0; count];
/// let lower_bound = -10.0;
/// let upper_bound = 10.0;
/// let step = (upper_bound - lower_bound) / (count - 1) as f64;
/// for i in 0..count {
///     let pos = lower_bound + i as f64 * step;
///     x[i] = pos;
///     y[i] = pos.sin();
/// }
/// let res = trapz_dt_sli(&x, &y);
/// assert!(res.abs() < 1.0e-3);
/// ```
pub fn trapz_dt_sli(x: &[f64], y: &[f64]) -> f64 {
    assert!(x.len() == y.len());
    let mut area = 0.0;
    for i in 0..x.len() - 1 {
        area += chunk_trapz(y[i], y[i + 1], x[i + 1] - x[i]);
    }
    area
}

/// # (d1) Simpson's rule - Data by slices
/// Integrate one-dimensional function represented by data enclosed in slices through the simpson rule:
/// * `x` - array of x-coordinates
/// * `y` - array of y-coordinates
/// 
/// ```
/// # use scilib::math::integration::d1::*;
/// let count = 11;
/// let mut x = vec![0.0; count];
/// let mut y = vec![0.0; count];
/// let lower_bound = -10.0;
/// let upper_bound = 10.0;
/// let step = (upper_bound - lower_bound) / (count - 1) as f64;
/// for i in 0..count {
///     let pos = lower_bound + i as f64 * step;
///     x[i] = pos;
///     y[i] = pos.sin();
/// }
/// let res = simp_dt_sli(&x, &y);
/// assert!(res.abs() < 1.0e-4);
/// ```
pub fn simp_dt_sli(x: &[f64], y: &[f64]) -> f64 {
    assert!(x.len() == y.len() && x.len() >= 3 && x.len() % 2 == 1);
    let mut area = 0.0;
    for i in (0..x.len() - 2).step_by(2) {
        area += chunk_simp_dt((x[i], y[i]), (x[i + 1], y[i + 1]), (x[i + 2], y[i + 2]));
    }
    area
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
