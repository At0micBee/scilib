//!
//! Implementation of the Bessel functions
//! 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

use std::f64::consts::{                 // Calling constants
    PI,
    FRAC_PI_2
};

use super::{                           // Calling sub-modules
    basic,
    complex::Complex
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # Precision limit for Bessel computation
const PRECISION_CONVERGENCE: f64 = 1.0e-8;

/// # Limit when computing Bessel Y
const DISTANCE_Y_LIM: f64 = 0.001;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # J Bessel function, integer index
///
/// `x` is the value to evaluate, and `n` the order of the function.
/// 
/// The J Bessel represent the first kind of Bessel function, and are solution to Bessel differential equation.
/// 
/// The bessel function depend on an infinite sum of terms; which we can't have.
/// The criterion chosen here is check each new term impacts the results significantly enough.
/// The default value selected in the program is defined by `const PRECISION_CONVERGENCE: f64 = 1.0e-8;`.
/// The J Bessel function cannot return values for `x < 0`, as they should be complex. We return 0 in those cases.
/// 
/// ```
/// # use scilib::math::complex::Complex;
/// # use scilib::math::bessel::j;
/// // Computing some example values
/// let res_00: Complex = j(0.0, 0);
/// let res_01: Complex = j(1.0, 0);
/// let res_10: Complex = j(0.0, 1);
/// let res_11: Complex = j(1.0, 1);
/// let res_20: Complex = j(0.0, 2);
/// let res_21: Complex = j(1.0, 2);
/// let res: Complex = j(5.2, 7);
/// 
/// // Comparing to tabulated data
/// assert_eq!(res_00, 1.0.into());
/// assert!((res_01.re - 0.7651976865).abs() < 1.0e-8);
/// assert_eq!(res_10, 0.0.into());
/// assert!((res_11.re - 0.44005058).abs() < 1.0e-8);
/// assert_eq!(res_20, 0.0.into());
/// assert!((res_21.re - 0.11490348).abs() < 1.0e-8);
/// assert!((res.re - 0.06544728).abs() < 1.0e-8);
/// 
/// // The function also handles negative orders
/// let pos1: Complex = j(3.2, 3);
/// let neg1: Complex = j(3.2, -3);
/// let pos2: Complex = j(2.45, 6);
/// let neg2: Complex = j(2.45, -6);
/// 
/// assert!(pos1 == -neg1);
/// assert!(pos2 == neg2);
/// 
/// // The input is treated as complex
/// let c: Complex = Complex::from(1, 2.5);
/// let res: Complex = j(c, 2);
/// ```
pub fn j<T: Into<Complex>>(x: T, n: i32) -> Complex {

    let np: i32 = n.abs();                                      // Getting the positive value of n
    let x2: Complex = x.into() / 2.0;                           // Halving x
    let mut k: i32 = 0;                                         // Order counter
    let mut d1: f64 = 1.0;                                      // First div
    let mut d2: f64 = basic::factorial(np as usize) as f64;     // Second div
    let mut sg: f64 = 1.0;                                      // Sign of the term

    let mut term: Complex = x2.powi(np) / d2;                   // The term at each step
    let mut res: Complex = Complex::default();                  // The result of the operation

    // If the first term is already too small we exit directly
    if term.modulus() < PRECISION_CONVERGENCE {
        return res;
    }

    // Computing the terms of the infinite series
    'convergence: loop {
        res += term;

        // If the changed compared to the final value is small we break
        if (term / res).modulus().abs() < PRECISION_CONVERGENCE {
            break 'convergence;
        }

        k += 1;                         // Incrementing value
        sg *= -1.0;                     // changing the sign of the term
        d1 *= k as f64;                 // Next value in the n! term
        d2 *= (np + k) as f64;          // Next value in the (n+k)! term
        term = sg * x2.powi(np + 2 * k) / (d1 * d2);
    }

    if n.is_negative() {
        (-1.0_f64).powi(np) * res
    } else {
        res
    }
}

/// # J Bessel function, real index
///
/// `x` is the value to evaluate, and `n` the order of the function.
/// 
/// Similar to the other J Bessel method, but this one allows the use of a real (float) index, rather
/// than an integer. This method is more costly to use than the other, and thus isn't recommended for
/// integer orders.
/// 
/// ```
/// # use scilib::math::complex::Complex;
/// # use scilib::math::bessel::{ j, jf };
/// // This method allows the computation of real index for J
/// let res_pos: Complex = jf(1.0, 2.5);
/// let res_neg: Complex = jf(2.4, -1.75);
/// assert!((res_pos.re - 0.04949681).abs() < 1.0e-4);
/// assert!((res_neg.re - 0.11990699).abs() < 1.0e-4);
/// 
/// // We can also check that the results are coherent with integers
/// let res_i: Complex = j(0.75, 2);
/// let res_r: Complex = jf(0.75, 2);
/// assert!((res_i.re - res_r.re).abs() < 1.0e-4);
/// 
/// // Because results are Complex, negative numbers are allowed
/// let neg: Complex = jf(-0.75, 2.3);
/// let expected: Complex = Complex::from(0.0219887007, 0.030264850);
/// 
/// assert!((neg.re - expected.re).abs() < 1.0e-5 && (neg.im - expected.im).abs() < 1.0e-5);
/// 
/// // As for j, we can also use Complex numbers
/// let c: Complex = Complex::from(1.2, 0.5);
/// let res: Complex = jf(c, 1.5);
/// 
/// assert!((res.re - 0.3124202913).abs() < 1.0e-5 && (res.im - 0.1578998151) < 1.0e-5);
/// ```
pub fn jf<T, U>(x: T, order: U) -> Complex
    where T: Into<Complex>, U: Into<f64> {

    let n: f64 = order.into();
    // If the number passed in whole, we fall back on the other method instead
    if n.fract() == 0.0 {
        return j(x, n as i32);
    }

    let x2: Complex = x.into() / 2.0;           // Halving x
    let mut k: f64 = 0.0;                       // Order counter
    let mut d1: f64 = 1.0;                      // First div
    let mut d2: f64 = basic::gamma(n + 1.0);    // Second div
    let mut sg: f64 = 1.0;                      // Sign of the term

    let mut term: Complex = x2.powf(n) / d2;    // The term at each step
    let mut res: Complex = Complex::default();  // The result of the operation
    
    // If the first term is already too small we exit directly
    if term.modulus().abs() < PRECISION_CONVERGENCE {
        return res;
    }

    // Computing the terms of the infinite series
    'convergence: loop {
        res += term;

        // If the changed compared to the final value is small we break
        if (term / res).modulus().abs() < PRECISION_CONVERGENCE {
            break 'convergence;
        }

        k += 1.0;                       // Incrementing value
        sg *= -1.0;                     // changing the sign of the term
        d1 *= k;                        // Next value in the n! term
        d2 *= n + k;                    // Next value in the gamma(n+k+1) term
        term = sg * x2.powf(n + 2.0 * k) / (d1 * d2);
    }

    res
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # Y Bessel function, real index
/// 
/// The Y Bessel represent the second kind of Bessel function, and are solution to Bessel differential equation,
/// in the case of a singularity at the origin.
/// 
/// `x` is the value to evaluate, and `n` the order of the function.
/// 
/// Because the function is not continuous for integer values of `n`, we need to compute the limit around these points.
/// We set the limit distance with `DISTANCE_Y_LIM`, compute the limit above and below the desired point and take the average.
/// We achieve precision under `1.0e-5` for non-integer`n`, and integer `n` using this approach. The functions are designed to
/// work with real numbers, thus the Y function cannot return a value for `x < 0`, as they become complex (returns 0).
/// 
/// ```
/// # use scilib::math::complex::Complex;
/// # use scilib::math::bessel::y;
/// let res_pos = y(1.0, 1.5);
/// let res_neg = y(1.0, -1.5);
/// 
/// assert!((res_pos.re - -1.10249557).abs() < 1.0e-5);
/// assert!((res_neg.re - -0.24029783).abs() < 1.0e-5);
/// 
/// // Values for integer n also works
/// let res_int_p = y(0.5, 1);
/// let res_int_n = y(0.5, -1);
/// 
/// assert!((res_int_p.re - -1.47147239).abs() < 1.0e-5);
/// assert!((res_int_n.re - 1.47147239).abs() < 1.0e-5);
/// 
/// // We can compute negative value with Y, the result is complex
/// let res_neg = y(-1.2, 3.1);
/// assert!((res_neg.re - 3.90596471).abs() < 1.0e-5 && (res_neg.im - -1.32157214).abs() < 1.0e-5);
/// 
/// // And we can use Complex as input
/// let c: Complex = Complex::from(-1.0, -0.5);
/// let res_c = y(c, 2.0);
/// 
/// assert!((res_c.re - -0.79108492).abs() < 1.0e-5 && (res_c.im - 0.60211151).abs() < 1.0e-5);
/// ```
pub fn y<T, U>(x: T, order: U) -> Complex
    where T: Into<Complex> + Copy, U: Into<f64> {

    let n: f64 = order.into();

    // If n is whole, we have to take the limit, otherwise it's direct
    if n.fract() == 0.0 {
        (y(x, n + DISTANCE_Y_LIM) + y(x, n - DISTANCE_Y_LIM)) / 2.0
    } else {
        ((n * PI).cos() * jf(x, n) - jf(x, -n)) / (n * PI).sin()
    }
}

/// # I modified Bessel function
pub fn i<T, U>(x: T, order: U) -> Complex
    where T: Into<Complex>, U: Into<f64> + Copy {
    
    // i^(-n) * jn(ix)
    Complex::i().powf(-order.into()) * jf(Complex::i() * x, order)
}

/// # K modified Bessel function
pub fn k<T, U>(x: T, order: U) -> Complex
    where T: Into<Complex> + Copy, U: Into<f64> {
    
    let n: f64 = order.into();

    let pos = i(x, n);
    let neg = i(x, -n);
    let factor: f64 = FRAC_PI_2 / (n * PI).sin();
    
    factor * (neg - pos)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # First Hankel function: H1
/// 
/// Computes the first kind of Hankel function, accepts complex input.
pub fn hankel_first<T, U>(x: T, order: U) -> Complex
    where T: Into<Complex> + Copy, U: Into<f64> {

    let n: f64 = order.into();
    let res_j = jf(x, n);
    let res_y = Complex::i() * y(x, n);

    res_j + res_y
}

/// # Second Hankel function: H2
/// 
/// Computes the second kind of Hankel function, accepts complex input.
/// We simplify the computation by simply conjugating the first kind
pub fn hankel_second<T, U>(x: T, order: U) -> Complex
    where T: Into<Complex> + Copy, U: Into<f64> {
    
    hankel_first(x, order).conjugate()
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
