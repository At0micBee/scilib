//!
//! # One-dimensional integration toolkit
//! 
//! One dimensional data integration functions.
//! 
//! ## Passing data to to the integrator
//! Integration can be done by directly passing slices with the x- and y-axis values, for example:
//! 
//! ```
//! # use scilib::integration::d1::*;
//! # use scilib::range;
//! let x = range::linear(0, 1, 1000);
//! let y: Vec<f64> = x.iter().map(|v| v.sin()).collect();
//! let res = trapeze(&x, &y);
//! assert!((res - 0.45969769413186028).abs() < 1.0e-7);
//! ```
//! 
//! Integrates the value of $\sin(x)$ between 0 and 1 with 1000 subdivisions using the trapeze method.
//! 
//! ## Passing the function as argument
//! Similarly, it is possible to integrate by directly passing the function (or closure):
//! 
//! ```
//! # use scilib::integration::d1::*;
//! let res = fn_simpson(|x| x.sin(), 0.0, 1.0, 1000);
//! assert!((res - 0.45969769413186028).abs() < 1.0e-7);
//! ```
//! 
//! Here, we also integrate $\sin(x)$ between 0 and 1 with 1000 subdivisions, but using the Simpson's method.
//! 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Internals

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
    
    (h0 + h1) / 6.0 * ((2.0 - h1 / h0) * p0.1 + (h0 + h1).powi(2) / (h0 * h1) * p1.1 + (2.0 - h0 / h1) * p2.1)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Passing function as argument

/// # Rectangle rule - Function as argument
/// 
/// ## Definition
/// Integration of a one-dimensional function using the left aligned rectangle rule:
/// $$
/// \int_a^b f(x)dx \approx \frac{b-a}{n}\sum_i^{n-1} f(x_i) 
/// $$
/// 
/// ## Inputs
/// - `function`: Closure or function pointer matching `f(x) = y`
/// - `lower_bound`: lower bound of the integral (Fixed value)
/// - `upper_bound`: Upper bound of the integral (Fixed value)
/// - `div`: Number of chunk evaluated
/// 
/// Returns the area under the curve
///
/// **CAUTION**: `div` must be greater than 1.
/// A _large_ value yields high precision & low performances,
/// while a _small_ value yields high performances & low precision
///
/// ## Example
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

/// # Trapezoidal rule - Function as argument
/// 
/// ## Definition
/// Integration of a one-dimensional function using the trapezoidal rule:
/// $$
/// \int_a^b f(x)dx \approx \frac{b-a}{n}\sum_i^{n-1} \frac{f(x_i) + f(x_{i+1})}{2}
/// $$
/// 
/// ## Inputs
/// - `function`: Closure or function pointer matching `f(x) = y`
/// - `lower_bound`: lower bound of the integral (Fixed value)
/// - `upper_bound`: Upper bound of the integral (Fixed value)
/// - `div`: Number of chunk evaluated
/// 
/// Returns the area under the curve
///
/// **CAUTION**: `div` must be greater than 1.
/// A _large_ value yields high precision & low performances,
/// while a _small_ value yields high performances & low precision
///
/// ## Example
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

/// # Simpson's rule - Function as argument
/// 
/// ## Definition
/// Integrate one-dimensional function through the simpson rule:
/// $$
/// \int_a^b f(x)dx \approx \frac{b-a}{3n}\sum_{i=1}^{n/2}
/// \left( f(x_{2i-2}) + 4f(x_{2i-1}) + f(x_{2i}) \right)
/// $$
/// 
/// ## Inputs
/// - `function`: Closure or function pointer matching `f(x) = y`
/// - `lower_bound`: lower bound of the integral (Fixed value)
/// - `upper_bound`: Upper bound of the integral (Fixed value)
/// - `div`: Number of chunk evaluated
/// 
/// Returns the area under the curve
///
/// **CAUTION**: `div` must be greater than 1.
/// A _large_ value yields high precision & low performances,
/// while a _small_ value yields high performances & low precision
///
/// ## Example
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

/// # Rectangle rule - X and Y Slices
/// 
/// ## Definition
/// Integration of a one-dimensional function using the left aligned rectangle rule:
/// $$
/// \int_a^b f(x)dx \approx \frac{b-a}{n}\sum_i^{n-1} f(x_i) 
/// $$
/// 
/// ## Inputs
/// Integrate one-dimensional function represented by data enclosed in slices through the rectangle rule:
/// - `x`: array of x-coordinates
/// - `y`: array of y-coordinates
/// 
/// Returns the area under the curve
///
/// **CAUTION**: `x.len()` and `y.len()` must be equal
///
/// ## Example
/// ```
/// # use scilib::integration::d1::*;
/// # use scilib::range;
/// let x = range::linear(0, 1, 1000);
/// let y: Vec<f64> = x.iter().map(|v| v.sin()).collect();
/// let res = rectangle(&x, &y);
/// assert!((res - 0.45969769413186028).abs() < 1.0e-3);
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

/// # Trapezoidal rule - X and Y Slices
/// 
/// ## Definition
/// Integration of a one-dimensional function using the trapezoidal rule:
/// $$
/// \int_a^b f(x)dx \approx \frac{b-a}{n}\sum_i^{n-1} \frac{f(x_i) + f(x_{i+1})}{2}
/// $$
/// 
/// ## Inputs
/// Integrate one-dimensional function represented by data enclosed in slices through the rectangle rule:
/// - `x`: array of x-coordinates
/// - `y`: array of y-coordinates
/// 
/// Returns the area under the curve
///
/// **CAUTION**: `x.len()` and `y.len()` must be equal
///
/// ## Example
/// ```
/// # use scilib::integration::d1::*;
/// # use scilib::range;
/// let x = range::linear(0, 1, 1000);
/// let y: Vec<f64> = x.iter().map(|v| v.sin()).collect();
/// let res = trapeze(&x, &y);
/// assert!((res - 0.45969769413186028).abs() < 1.0e-7);
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

/// # Simpson's rule - X and Y Slices
/// 
/// ## Definition
/// Integrate one-dimensional function through the simpson rule:
/// $$
/// \int_a^b f(x)dx \approx \frac{b-a}{3n}\sum_{i=1}^{n/2}
/// \left( f(x_{2i-2}) + 4f(x_{2i-1}) + f(x_{2i}) \right)
/// $$
/// 
/// ## Inputs
/// Integrate one-dimensional function represented by data enclosed in slices through the rectangle rule:
/// - `x`: array of x-coordinates
/// - `y`: array of y-coordinates
/// 
/// Returns the area under the curve
///
/// **CAUTION**: `x.len()` and `y.len()` must be equal
///
/// ## Example
/// ```
/// # use scilib::integration::d1::*;
/// # use scilib::range;
/// let x = range::linear(0, 1, 1001);
/// let y: Vec<f64> = x.iter().map(|v| v.sin()).collect();
/// let res = trapeze(&x, &y);
/// assert!((res - 0.45969769413186028).abs() < 1.0e-7);
/// ```
pub fn simpson(x: &[f64], y: &[f64]) -> f64 {

    assert!(
        x.len() == y.len() && x.len() >= 3 && x.len() % 2 == 1,
        "`x.len()` and `y.len()` must be equals, and be an odd number greater than 3!"
    );

    (0..(x.len() - 1)).step_by(2).fold(0.0, |sum, idx| {
        sum + chunk_simp_dt((x[idx], y[idx]), (x[idx + 1], y[idx + 1]), (x[idx + 2], y[idx + 2]))
    })
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// X and Y as tuples

/// # Rectangle rule - Data tuple
/// 
/// ## Definition
/// Integration of a one-dimensional function using the left aligned rectangle rule:
/// $$
/// \int_a^b f(x)dx \approx \frac{b-a}{n}\sum_i^{n-1} f(x_i) 
/// $$
/// 
/// ## Inputs
/// - `data`: array of tuples corresponding to `x` and `y`, respectively
///
/// Returns the area under the curve
/// 
/// ## Example
/// ```
/// # use scilib::integration::d1::*;
/// # use scilib::range;
/// let x = range::linear(0, 1, 1000);
/// let y: Vec<f64> = x.iter().map(|v| v.sin()).collect();
/// let data: Vec<(f64, f64)> = x.into_iter().zip(y.into_iter()).map(|(a, b)| {
///     (a, b)
/// }).collect();
/// let res = tuple_rectangle(&data);
/// assert!((res - 0.45969769413186028).abs() < 1.0e-3);
/// ```
pub fn tuple_rectangle(data: &[(f64, f64)]) -> f64 {

    (0..(data.len() - 1)).fold(0.0, |sum, idx| {
        sum + chunk_rect(data[idx].1, data[idx + 1].0 - data[idx].0)
    })
}

/// # Trapezoidal rule - Data tuple
/// 
/// ## Definition
/// Integration of a one-dimensional function using the trapezoidal rule:
/// $$
/// \int_a^b f(x)dx \approx \frac{b-a}{n}\sum_i^{n-1} \frac{f(x_i) + f(x_{i+1})}{2}
/// $$
/// 
/// ## Inputs
/// - `data`: array of tuples corresponding to `x` and `y`, respectively
///
/// Returns the area under the curve
/// 
/// ## Example
/// ```
/// # use scilib::integration::d1::*;
/// # use scilib::range;
/// let x = range::linear(0, 1, 1000);
/// let y: Vec<f64> = x.iter().map(|v| v.sin()).collect();
/// let data: Vec<(f64, f64)> = x.into_iter().zip(y.into_iter()).map(|(a, b)| {
///     (a, b)
/// }).collect();
/// let res = tuple_trapeze(&data);
/// assert!((res - 0.45969769413186028).abs() < 1.0e-7);
/// ```
pub fn tuple_trapeze(data: &[(f64, f64)]) -> f64 {

    (0..(data.len() - 1)).fold(0.0, |sum, idx| {
        sum + chunk_trapz(data[idx].1, data[idx + 1].1, data[idx + 1].0 - data[idx].0)
    })
}

/// # Simpson's rule - Data tuple
/// 
/// ## Definition
/// Integrate one-dimensional function through the simpson rule:
/// $$
/// \int_a^b f(x)dx \approx \frac{b-a}{3n}\sum_{i=1}^{n/2}
/// \left( f(x_{2i-2}) + 4f(x_{2i-1}) + f(x_{2i}) \right)
/// $$
/// 
/// ## Inputs
/// - `data`: array of tuples corresponding to `x` and `y`, respectively
///
/// Returns the area under the curve
/// 
/// ## Example
/// ```
/// # use scilib::integration::d1::*;
/// # use scilib::range;
/// let x = range::linear(0, 1, 1001);
/// let y: Vec<f64> = x.iter().map(|v| v.sin()).collect();
/// let data: Vec<(f64, f64)> = x.into_iter().zip(y.into_iter()).map(|(a, b)| {
///     (a, b)
/// }).collect();
/// let res = tuple_simpson(&data);
/// assert!((res - 0.45969769413186028).abs() < 1.0e-7);
/// ```
pub fn tuple_simpson(data: &[(f64, f64)]) -> f64 {

    assert!(
        data.len() >= 3 && data.len() % 2 == 1,
        "size of data parameter must be an uneven number greater than 3!"
    );

    (0..(data.len() - 1)).step_by(2).fold(0.0, |sum, idx| {
        sum + chunk_simp_dt(data[idx], data[idx + 1], data[idx + 2])
    })
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
