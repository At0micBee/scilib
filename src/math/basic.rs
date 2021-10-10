//!
//! # Basic math functions
//! 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

use super::super::constant;     // Calling the constants from the module

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// Precision used for convergence
const PRECISION: f64 = 1.0e-12;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # Sinus cardinal
///
/// `x`: the value at which to evaluate the function
/// 
/// Returns: the sinc value of x
/// 
/// ```
/// # use abacus::math::basic::sinc;
/// let x: f64 = sinc(0.0);
/// let y: f64 = sinc(1.2);
/// 
/// // Comparing to tabulated values
/// assert_eq!(x, 1.0);
/// assert!((y - 0.776699238306) < 1.0e-12);
/// ```
pub fn sinc(x: f64) -> f64 {
    let res = x.sin() / x;
    if res.is_nan() {
        1.0
    } else {
        res
    }
}

/// # Newton binomial formula
///
/// `n` is the number of "objects" and `k` is the selection
/// 
/// Returns: `k` among `n`
/// 
/// ```
/// # use abacus::math::basic::binomial;
/// let res: usize = binomial(4_usize, 2_usize);
/// assert_eq!(res, 6);
/// ```
/// 
/// The implementation used here relies on the recurrence relation:
/// 
/// (n, k) = (n-1, k) + (n-1, k-1)
/// 
/// This pushes overflow back for a few more terms.
pub fn binomial(n: usize, k: usize) -> usize {

    // n must be greater than k to produce a value
    if k > n {
        return 0;
    }

    let mut res: usize = 1;             // We initialize the result at 1
    let mut n_d: usize = n;             // We need a mutable value for n

    // We loop the counter up to k possible values
    for val in 1..=k {
        res *= n_d;     // Changing the result by n
        n_d -= 1;       // Decrementing n
        res /= val;     // Changing the result by the choices options
    }

    // Returning the result
    res
}

/// # Factorial function
///
/// `n` is the integer at which to evaluate the factorial
/// 
/// Returns: `n`!
/// 
/// ```
/// # use abacus::math::basic::factorial;
/// let res: usize = factorial(5_usize);
/// assert_eq!(res, 120);
/// ```
pub fn factorial<T>(n: T) -> usize
    where T: Into<usize> {
    (1..=n.into()).fold(1, |res, val| res * val)
}

/// # Gamma function
/// 
/// `x` is the value to evaluate.
/// 
/// The gamma function is the generalization of the factorial function, and is
/// used to provide value for non-integer numbers (except non-positive integers).
/// 
/// ```
/// # use abacus::math::basic::gamma;
/// let res_1: f64 = gamma(2.3);
/// let res_2: f64 = gamma(-0.45);
/// 
/// assert!((res_1 - 1.16671190).abs() < 1.0e-5);
/// assert!((res_2 - -3.591387).abs() < 1.0e-5);
/// ```
/// 
/// With the current computation scheme, we limit the precision of the computation in exchange for speed.
/// Typical values are achieve within a `1.0e-5` margin of error. Changing the method to another one
/// might grant some more speed and lower the error on the results.
pub fn gamma(x: f64) -> f64 {

    let mut n: f64 = 1.0;      // Order counter

    // The values of each term and the result
    let mut term: f64 = x.exp() / (1.0 + x);
    let mut res: f64 = 1.0;

    // If the first term is already too small we exit directly
    if (term - 1.0).abs() < PRECISION {
        return res;
    }

    // Computing the terms of the infinite series
    'convergence: loop {
        res *= term;

        //If the changes become too small, we stop
        if (term - 1.0).abs() < PRECISION {
            break 'convergence;
        }

        // Updating the values
        n += 1.0;
        term = (x / n).exp() / (1.0 + x / n);
    }

    res * (-x * constant::EULER_MASCHERONI).exp() / x
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
