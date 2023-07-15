//!
//! # One-dimensional integration toolkit
//! 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[inline(always)]
fn chunk_rect(a: f64, offset: f64) -> f64 {
    a * offset
}

#[inline(always)]
fn chunk_trapz(a: f64, b: f64, offset: f64) -> f64 {
    (a + b) / 2.0 * offset
}

#[inline(always)]
fn chunk_simp(
        p0: (f64 /*a0*/, f64 /*b0*/),
        p1: f64, /*b1*/
        p2: (f64 /*a2*/, f64 /*b2*/),
    ) -> f64 {

    (p2.0 - p0.0) / 6.0 * (p0.1 + 4.0 * p1 + p2.1)
}

#[inline(always)]
fn chunk_simp_dt(
        p0: (f64 /*a0*/, f64 /*b0*/),
        p1: (f64 /*a1*/, f64 /*b1*/),
        p2: (f64 /*a2*/, f64 /*b2*/),
    ) -> f64 {

    let h0 = p1.0 - p0.0;
    let h1 = p2.0 - p1.0;
    
    (h0 + h1) / 6.0
        * ((2.0 - h1 / h0) * p0.1 + (h0 + h1).powi(2) / (h0 * h1) * p1.1 + (2.0 - h0 / h1) * p2.1)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Passing function as argument

/// # (d1) Rectangle rule - Function
/// Integrate one-dimensional function through the rectangle rule:
/// * `function` - Closure or function pointer matching `f(x) = y`
/// * `lower_bound` - First limit of the integral (Fixed value)
/// * `upper_bound` - Second limit of the integral (Fixed value)
/// * `div` - Number of chunk evaluated: _big_ = high precision & low performances, _small_ = high performances & low precision
///
/// **CAUTION**: `div` must be greater than 1
///
/// ```
/// # use scilib::integration::d1::*;
/// # use std::f64::consts::PI;
/// let res = fn_rectangle(|x| x.cos(), PI * 2.0, PI * 2.0 + PI / 2.0, 1000);
/// assert!((res - 1.0).abs() < 1.0e-3);
/// ```
pub fn fn_rectangle(
        function: impl Fn(f64) -> f64,
        lower_bound: f64,
        upper_bound: f64,
        div: usize,
    ) -> f64 {

    assert!(div >= 1, "div parameter must be greater than 1!");

    let step = (upper_bound - lower_bound) / div as f64;
    
    (0..div).fold(0.0, |sum, idx| {
        sum + chunk_rect(function(step * (idx as f64) + lower_bound), step)
    })
}

/// # (d1) Trapezoidal rule - Function
/// Integrate one-dimensional function through the trapezoidal rule:
/// * `function` - Closure or function pointer matching `f(x) = y`
/// * `lower_bound` - First limit of the integral (Fixed value)
/// * `upper_bound` - Second limit of the integral (Fixed value)
/// * `div` - Number of chunk evaluated: _big_ = high precision & low performances, _small_ = high performances & low precision
///
/// **CAUTION**: `div` must be greater than 1
///
/// ```
/// # use scilib::integration::d1::*;
/// let res = fn_trapeze(|x| x.exp(), -2.0, 10.0, 10_000);
/// assert!((res - 22026.3304).abs() < 1.0e-2);
/// ```
pub fn fn_trapeze(
        function: impl Fn(f64) -> f64,
        lower_bound: f64,
        upper_bound: f64,
        div: usize,
    ) -> f64 {

    assert!(div >= 1, "div parameter must be greater than 1!");

    let step = (upper_bound - lower_bound) / div as f64;
    let mut lower: f64 = function(lower_bound);
    let mut upper: f64 = lower;

    (0..div).fold(0.0, |mut sum, idx| {
        upper = function(step * (idx as f64 + 1.0) + lower_bound);
        sum += chunk_trapz(lower, upper, step);
        lower = upper;
        sum
    })
}

/// # (d1) Simpson's rule - Function
/// Integrate one-dimensional function through the simpson rule:
/// * `function` - Closure or function pointer matching `f(x) = y`
/// * `lower_bound` - First limit of the integral (Fixed value)
/// * `upper_bound` - Second limit of the integral (Fixed value)
/// * `div` - Number of chunk evaluated: _big_ = high precision & low performances, _small_ = high performances & low precision
///
/// **CAUTION**: `div` must be an even number greater than 2
///
/// ```
/// # use scilib::integration::d1::*;
/// let res = fn_simpson(|x| x.exp(), -2.0, 10.0, 1000);
/// assert!((res - 22026.330459).abs() < 1.0e-4);
/// ```
pub fn fn_simpson(
        function: impl Fn(f64) -> f64,
        lower_bound: f64,
        upper_bound: f64,
        div: usize,
    ) -> f64 {

    assert!(
        div >= 2 && div % 2 == 0,
        "div parameter must be an even number greater than 2!"
    );

    let step = (upper_bound - lower_bound) / div as f64;
    let mut p0: (f64, f64) = (lower_bound, function(lower_bound));
    let mut p2: (f64, f64) = (0.0, 0.0);
    let mut curr_pos: f64 = 0.0;

    (0..div).step_by(2).fold(0.0, |mut sum, idx| {
        curr_pos = lower_bound + (idx + 2) as f64 * step;
        p2 = (curr_pos, function(curr_pos));
        sum += chunk_simp(p0, function((p0.0 + p2.0) / 2.0), p2);
        p0 = p2;
        sum
    })
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// X and Y vectors as input

/// # (d1) Rectangle rule - Data by slices
/// Integrate one-dimensional function represented by data enclosed in slices through the rectangle rule:
/// * `x` - array of x-coordinates
/// * `y` - array of y-coordinates
///
/// **CAUTION**: `x.len()` and `y.len()` must be equals
///
/// ```
/// # use scilib::integration::d1::*;
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
/// let res = rectangle(&x, &y);
/// assert!(res.abs() < 1.0e-1);
/// ```
pub fn rectangle(x: &[f64], y: &[f64]) -> f64 {

    assert!(
        x.len() == y.len(),
        "`x.len()` and `y.len()` must be equals!"
    );

    (0..(x.len() - 1)).fold(0.0, |sum, idx| {
        sum + chunk_rect(y[idx], x[idx + 1] - x[idx])
    })
}

/// # (d1) Trapezoidal rule - Data by slices
/// Integrate one-dimensional function represented by data enclosed in slices through the trapezoidal rule:
/// * `x` - array of x-coordinates
/// * `y` - array of y-coordinates
///
/// **CAUTION**: `x.len()` and `y.len()` must be equals
///
/// ```
/// # use scilib::integration::d1::*;
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
/// let res = trapeze(&x, &y);
/// assert!(res.abs() < 1.0e-3);
/// ```
pub fn trapeze(x: &[f64], y: &[f64]) -> f64 {

    assert!(
        x.len() == y.len(),
        "`x.len()` and `y.len()` must be equals!"
    );

    (0..(x.len() - 1)).fold(0.0, |sum, idx| {
        sum + chunk_trapz(y[idx], y[idx + 1], x[idx + 1] - x[idx])
    })
}

/// # (d1) Simpson's rule - Data by slices
/// Integrate one-dimensional function represented by data enclosed in slices through the simpson rule:
/// * `x` - array of x-coordinates
/// * `y` - array of y-coordinates
///
/// **CAUTION**: `x.len()` and `y.len()` must be equals and `x.len()` must be an uneven number greater than 3
///
/// ```
/// # use scilib::integration::d1::*;
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
/// let res = simpson(&x, &y);
/// assert!(res.abs() < 1.0e-4);
/// ```
pub fn simpson(x: &[f64], y: &[f64]) -> f64 {

    assert!(
        x.len() == y.len() && x.len() >= 3 && x.len() % 2 == 1,
        "`x.len()` and `y.len()` must be equals and `x.len()` must be an uneven number greater than 3!"
    );

    (0..(x.len() - 1)).step_by(2).fold(0.0, |sum, idx| {
        sum + chunk_simp_dt((x[idx], y[idx]), (x[idx + 1], y[idx + 1]), (x[idx + 2], y[idx + 2]))
    })
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// X and Y as tuples

/// # (d1) Rectangle rule - Data by tuples
/// Integrate one-dimensional function represented by data enclosed in tuples through the rectangle rule:
/// * `data` - array of tuples corresponding to `x` and `y` respectively
///
/// ```
/// # use scilib::integration::d1::*;
/// let count = 1_000;
/// let mut data = vec![(0.0, 0.0); count];
/// let lower_bound = -10.0;
/// let upper_bound = 10.0;
/// let step = (upper_bound - lower_bound) / (count - 1) as f64;
/// for i in 0..count {
///     let x = lower_bound + i as f64 * step;
///     data[i] = (x, x.sin());
/// }
/// let res = tuple_rectangle(&data);
/// assert!(res.abs() < 1.0e-1);
/// ```
pub fn tuple_rectangle(data: &[(f64, f64)]) -> f64 {

    (0..(data.len() - 1)).fold(0.0, |sum, idx| {
        sum + chunk_rect(data[idx].1, data[idx + 1].0 - data[idx].0)
    })
}

/// # (d1) Trapezoidal rule - Data by tuples
/// Integrate one-dimensional function represented by data enclosed in tuples through the trapezoidal rule:
/// * `data` - array of tuples corresponding to `x` and `y` respectively
///
/// ```
/// # use scilib::integration::d1::*;
/// let count = 100;
/// let mut data = vec![(0.0, 0.0); count];
/// let lower_bound = -10.0;
/// let upper_bound = 10.0;
/// let step = (upper_bound - lower_bound) / (count - 1) as f64;
/// for i in 0..count {
///     let x = lower_bound + i as f64 * step;
///     data[i] = (x, x.sin());
/// }
/// let res = tuple_trapeze(&data);
/// assert!(res.abs() < 1.0e-3);
/// ```
pub fn tuple_trapeze(data: &[(f64, f64)]) -> f64 {

    (0..(data.len() - 1)).fold(0.0, |sum, idx| {
        sum + chunk_trapz(data[idx].1, data[idx + 1].1, data[idx + 1].0 - data[idx].0)
    })
}

/// # (d1) Simpson's rule - Data by tuples
/// Integrate one-dimensional function represented by data enclosed in tuples through the trapezoidal rule:
/// * `data` - array of tuples corresponding to `x` and `y` respectively
///
/// **CAUTION**: `data.len()` must be an uneven number greater than 3
///
/// ```
/// # use scilib::integration::d1::*;
/// let count = 11;
/// let mut data = vec![(0.0, 0.0); count];
/// let lower_bound = -10.0;
/// let upper_bound = 10.0;
/// let step = (upper_bound - lower_bound) / (count - 1) as f64;
/// for i in 0..count {
///     let x = lower_bound + i as f64 * step;
///     data[i] = (x, x.sin());
/// }
/// let res = tuple_simpson(&data);
/// assert!(res.abs() < 1.0e-4);
/// ```
pub fn tuple_simpson(data: &[(f64, f64)]) -> f64 {

    assert!(
        data.len() >= 3 && data.len() % 2 == 1,
        "size of data paramater must be an uneven number greater than 3!"
    );

    (0..(data.len() - 1)).step_by(2).fold(0.0, |sum, idx| {
        sum + chunk_simp_dt(data[idx], data[idx + 1], data[idx + 2])
    })
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
