//!
//! # Bessel functions
//! 
//! The Bessel functions are the solution to the [Bessel differential equation](https://en.wikipedia.org/wiki/Bessel_function#Spherical_Bessel_functions).
//! There are multiple variants of these solutions, and this sub-module provides functions for all of them.
//! 
//! ## First kind: J
//! 
//! The Jn(x) function is the simplest solution to Bessel's equation. The most generalized case is for
//! a real index `n` and a complex number `x`, which is possible with this crate. Additionally, two methods are
//! provided for this equation, one for a integer order (named `j`) and one for a real order (named `jf`). The integer
//! order one is limited but is faster than it's counterpart. The real order function is slower, but can compute any
//! order. Note that `jf` will attempt to fall back on `j` when it finds an integer order.
//! 
//! ```rust
//! # use scilib::math::complex::Complex;
//! # use scilib::math::bessel::{ j, jf };
//! let c = Complex::from(-0.75, 3);
//! let res_i = j(c, 1);                // Faster for integer order
//! let res_f = jf(c, 1.5);             // Would also work with 1.0
//! ```
//! 
//! ## Second kind: Y
//! 
//! Similar to the first kind, the Y equation are solution of Bessel's equation with a singularity at the origin.
//! The Y function is itself based on the J function. The function is undefined for any integer order, in which
//! case the limit has to be taken.
//! 
//! ```rust
//! # use scilib::math::complex::Complex;
//! # use scilib::math::bessel::y;
//! let c = Complex::from(2, -1.2);
//! let res_f = y(c, 1.5);              // Not a problem
//! let res_i = y(c, 1);                // The function takes the limit in this case
//! ```
//! 
//! ## Modified first kind: I
//! 
//! Also known as the hyperbolic Bessel function of the first kind. Its definition is similar to J, but lacks the
//! alternating `(-1)^k` term in the sum.
//! 
//! ```rust
//! # use scilib::math::complex::Complex;
//! # use scilib::math::bessel::i;
//! let c = Complex::from(0.2, 1);
//! let res = i(c, -1.2);
//! ```
//! 
//! ## Modified second kind: K
//! 
//! Also known as the hyperbolic Bessel function of the second kind. Its definition is similar to Y, but lacks the
//! `cos(n*pi)`, and is normalized by `pi/2`.
//! 
//! ```rust
//! # use scilib::math::complex::Complex;
//! # use scilib::math::bessel::k;
//! let c = Complex::from(0, 7);
//! let res = k(c, 0);
//! ```
//! 
//! ## Hankel functions: H1 and H2
//! 
//! Hankel functions are two linearly independent solutions to Bessel's equation.
//! 
//! ```rust
//! # use scilib::math::complex::Complex;
//! # use scilib::math::bessel::{ hankel_first, hankel_second };
//! let c = Complex::from(-0.3, 1.52);
//! let res_1 = hankel_first(c, -2.3);
//! let res_2 = hankel_second(c, -2.3);
//! ```
//! 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

use std::f64::consts::{     // Using std lib constants
    PI,                     // Pi
    FRAC_PI_2               // Pi / 2
};

use super::{                // Using parts from the crate
    basic,                  // Basic functions
    complex::Complex        // Using Complex numbers
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
/// integer orders. The function tries to prevent this by trying to fall back on the integer order version
/// when possible, but could fail in edge cases.
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
/// We achieve precision under `1.0e-5` for non-integer`n`, and integer `n` using this approach.
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
/// 
/// The I modified First Bessel function represent another kind of solution to the Bessel differential equation.
/// 
/// `x` is the value to evaluate, and `n` the order of the function.
/// 
/// We use a definition of I based on an infinite series (similar to J). This way, we ensure good precision in
/// the computation.
/// 
/// ```
/// # use scilib::math::complex::Complex;
/// # use scilib::math::bessel::i;
/// let res = i(1.2, 0);
/// assert!((res.re - 1.39373).abs() < 1.0e-4 && res.im == 0.0);
/// 
/// let c = Complex::from(-1.2, 0.5);
/// let r2 = i(c, -1.6);
/// assert!((r2.re - 0.549831).abs() < 1.0e-5 && (r2.im - -0.123202).abs() < 1.0e-5);
/// ```
pub fn i<T, U>(x: T, order: U) -> Complex
where T: Into<Complex>, U: Into<f64> + Copy {
    
    let n: f64 = order.into();

    let x2: Complex = x.into() / 2.0;           // Halving x
    let mut k: f64 = 0.0;                       // Order counter
    let mut d1: f64 = 1.0;                      // First div
    let mut d2: f64 = basic::gamma(n + 1.0);    // Second div

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
        d1 *= k;                        // Next value in the n! term
        d2 *= n + k;                    // Next value in the gamma(n+k+1) term
        term = x2.powf(n + 2.0 * k) / (d1 * d2);
    }

    res
}

/// # K modified Bessel function
/// 
/// The K modified Second Bessel function represent another kind of solution to the Bessel differential equation.
/// 
/// `x` is the value to evaluate, and `n` the order of the function.
/// 
/// The definition of K is similar to Y, but is based on I and not J.
/// 
/// ```
/// # use scilib::math::complex::Complex;
/// # use scilib::math::bessel::k;
/// let c1 = Complex::from(2, -1);
/// let res = k(c1, -3.5);
/// assert!((res.re - -0.32113627).abs() < 1.0e-5 && (res.im - 0.76751785).abs() < 1.0e-5);
/// 
/// // Similar to Y, we take the limit for integer orders
/// let c2 = Complex::from(-1.1, 0.6);
/// let res_i = k(c2, 1);
/// assert!((res_i.re - -1.6153940).abs() < 1.0e-5 && (res_i.im - -2.1056846).abs() < 1.0e-5);
/// ```
pub fn k<T, U>(x: T, order: U) -> Complex
where T: Into<Complex> + Copy, U: Into<f64> {

    let n: f64 = order.into();

    // If n is whole, we have to take the limit, otherwise it's direct
    if n.fract() == 0.0 {
        (k(x, n + DISTANCE_Y_LIM) + k(x, n - DISTANCE_Y_LIM)) / 2.0
    } else {
        (FRAC_PI_2 / (n * PI).sin()) * (i(x, -n) - i(x, n))
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # First Hankel function: H1
/// 
/// Computes the first kind of Hankel function, accepts complex input.
/// 
/// ```
/// # use scilib::math::complex::Complex;
/// # use scilib::math::bessel::hankel_first;
/// let c1 = Complex::from(-1.1, 2.3);
/// let r1 = hankel_first(c1, 1);
/// assert!((r1.re - -0.0112027).abs() < 1.0e-5 && (r1.im - 0.0551947).abs() < 1.0e-5);
/// 
/// let c2 = Complex::from(5.2, -3);
/// let r2 = hankel_first(c2, -2.35);
/// assert!((r2.re - -4.2809477).abs() < 1.0e-5 && (r2.im - 3.2123502).abs() < 1.0e-5);
/// ```
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
/// 
/// ```
/// # use scilib::math::complex::Complex;
/// # use scilib::math::bessel::hankel_second;
/// let c1 = Complex::from(-1.1, 2.3);
/// let r1 = hankel_second(c1, 1);
/// assert!((r1.re - -3.54421).abs() < 1.0e-5 && (r1.im - 2.2983539).abs() < 1.0e-5);
/// 
/// let c2 = Complex::from(5.2, -3);
/// let r2 = hankel_second(c2, -2.35);
/// assert!((r2.re - -0.0068184520).abs() < 1.0e-5 && (r2.im - -0.0193698).abs() < 1.0e-5);
pub fn hankel_second<T, U>(x: T, order: U) -> Complex
where T: Into<Complex> + Copy, U: Into<f64> {
    
    let n: f64 = order.into();
    let res_j = jf(x, n);
    let res_y = Complex::i() * y(x, n);

    res_j - res_y
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # First spherical Bessel function: j
/// 
/// Compute the first kind of spherical bessel function.
/// * `z` - where the function is evaluated
/// * `n` - order evaluated
/// 
/// ```
/// # use scilib::math::complex::Complex;
/// # use scilib::math::bessel::*;
/// 
/// let res = sj(Complex::from(13, 5), 3);
/// assert_eq!(res, Complex::from(1.6109825049200244, -4.322382277910797));
/// ```
pub fn sj<T>(z: T, n: usize) -> Complex 
where T: Into<Complex> {
    let x: Complex  = z.into();
    (std::f64::consts::PI / 2.0 / x).sqrt() * jf(x, n as f64 + 0.5)
}

fn sj_upward_reccurence<T>(z: T, n: usize) -> Vec<Complex> 
where T: Into<Complex> {
    let count = n + 1;
    let x: Complex = z.into();
    let mut jn = vec![Complex::from(0.0, 0.0); count];
    jn[0] = x.sin() / x;
    jn[1] = x.sin() / x.powi(2) - x.cos() / x;
    for i in 1..=count - 2 {
        jn[i + 1] = (2 * i + 1) as f64 / x * jn[i] - jn[i - 1];
    }
    jn
}

fn sj_downward_reccurence<T>(z: T, nl: usize, nu: usize) -> Vec<Complex> 
where T: Into<Complex> {
    let count = nu - nl + 1;
    let x: Complex = z.into();
    let mut jn= vec![Complex::from(0.0, 0.0); count + 2];
    jn[count + 1] = Complex::from(0.0, 0.0);
    jn[count] = Complex::from(1.0, 0.0);
    for i in (1..=count).rev() {
        jn[i - 1] = (2 * i + 1) as f64 / x * jn[i] - jn[i + 1];
    }
    jn.resize(count, Complex::from(0.0, 0.0));
    jn
}

/// # First spherical Bessel function (array): j
/// 
/// Compute the first kind of spherical bessel function by reccurency to get an array of function.
/// * `z` - where the function is evaluated
/// * `n` - maximum order evaluated, return a vector from 0 to the nth order included : `sj[0]` to `sj[n]`
/// 
/// ```
/// # use scilib::math::complex::Complex;
/// # use scilib::math::bessel::*;
/// 
/// let res1 = sj_array(Complex::from(13, 5), 3);
/// assert_eq!(res1[0], Complex::from(3.8248700377925635, 3.708547263317134));
/// assert_eq!(res1[1], Complex::from(-3.357143112679857, 3.9747696875545517));
/// assert_eq!(res1[2], Complex::from(-4.1924320794482135, -2.6499227040139104));
/// let res2 = sj_array(0.2, 25);
/// assert_eq!(res2[13], Complex::from(3.835110596379198e-24, 0.0));
/// assert_eq!(res2[17], Complex::from(5.910455642760406e-33, 0.0));
/// assert_eq!(res2[25], Complex::from(1.125476749298975e-51, 0.0));
/// ```
pub fn sj_array<T>(z: T, n: usize) -> Vec<Complex>
where T: Into<Complex> {
    let x: Complex = z.into();
    if x.modulus() > n as f64 / 2.0 {
        sj_upward_reccurence(x, n)
    } else {
        const PACK: usize = 50;
        let num_big_loop = n / PACK;
        let mut jn_all = Vec::<Vec<Complex>>::new();
        for i in 0..num_big_loop {
            jn_all.push(sj_downward_reccurence(x, i * PACK, (i + 1) * PACK));
        }
        let rest = n % PACK;
        if rest != 0 {
            jn_all.push(sj_downward_reccurence(x, n - rest, n));
        }

        let mut jn = Vec::<Complex>::with_capacity(n);
        let mut norm = x.sin() / x / jn_all[0][0];
        for i in 0..jn_all[0].len() {
            jn.push(jn_all[0][i] * norm);
        }
        for i in 1..jn_all.len() {
            norm = *jn.last().unwrap() / jn_all[i][0];
            for k in 1..jn_all[i].len() {
                jn.push(jn_all[i][k] * norm);
            }
        }
        jn
    }
}

/// # Second spherical Bessel function: y
/// 
/// Compute the second kind of spherical bessel function.
/// * `z` - where the function is evaluated
/// * `n` - order evaluated
/// 
/// ```
/// # use scilib::math::complex::Complex;
/// # use scilib::math::bessel::*;
/// 
/// let res = sy(Complex::from(13, 5), 3);
/// assert_eq!(res, Complex::from(4.322629120777188, 1.6104674841102558));
/// ```
pub fn sy<T>(z: T, n: usize) -> Complex 
where T: Into<Complex> {
    let x: Complex = z.into();
    (std::f64::consts::PI / 2.0 / x).sqrt() * y(x, n as f64 + 0.5)
}

/// # Second spherical Bessel function (array): y
/// 
/// Compute the second kind of spherical bessel function by reccurency to get an array of function.
/// * `z` - where the function is evaluated
/// * `n` - maximum order evaluated, return a vector from 0 to the nth order included : `sy[0]` to `sy[n]`
/// 
/// ```
/// # use scilib::math::complex::Complex;
/// # use scilib::math::bessel::*;
/// 
/// let res1 = sy_array(Complex::from(13, 5), 3);
/// assert_eq!(res1[0], Complex::from(-3.7090299518957797, 3.8248379131516654));
/// assert_eq!(res1[1], Complex::from(-3.9748349852610523, -3.356650136356049));
/// assert_eq!(res1[2], Complex::from(2.6504303824601, -4.1922958025278));
/// let res2 = sy_array(0.2, 25);
/// assert_eq!(res2[13], Complex::from(-4.82921204481494e22, 0.0));
/// assert_eq!(res2[17], Complex::from(-2.417182573235861e31, 0.0));
/// assert_eq!(res2[25], Complex::from(-8.711173815326792e49, 0.0));
/// ```
pub fn sy_array<T>(z: T, n: usize) -> Vec<Complex>
where T: Into<Complex> {
    let x: Complex = z.into();
    let y0 = -x.cos() / x;
    if n == 0 {
        return vec![y0];
    }
    let y1 = -x.cos() / (x * x) - x.sin() / x;
    if n == 1 {
        return vec![y0, y1];
    }
    let mut yn: Vec<Complex> = vec![Complex::from(0.0, 0.0); n + 1];
    yn[0] = y0;
    yn[1] = y1;
    for i in 1..=n - 1 {
        yn[i + 1] = ((2 * i + 1) as f64) / x * yn[i] - yn[i - 1];
    }
    return yn;
}

/// # First spherical Hankel function: h1
/// 
/// Compute the first kind of spherical hankel function.
/// * `z` - where the function is evaluated
/// * `n` - order evaluated
/// 
/// ```
/// # use scilib::math::complex::Complex;
/// # use scilib::math::bessel::*;
/// 
/// let res = sh_first(Complex::from(13, 5), 3);
/// assert_eq!(res, Complex::from(5.150208097686182e-4, 2.4684286639153896e-4));
/// ```
pub fn sh_first<T>(z: T, n: usize) -> Complex 
where T: Into<Complex> {
    let x: Complex = z.into();
    sj(x, n) + Complex::i() * sy(x, n)
}

/// # First spherical Hankel function (array): h1
/// 
/// Compute the first kind of spherical hankel function by reccurency to get an array of function.
/// * `z` - where the function is evaluated
/// * `n` - maximum order evaluated, return a vector from 0 to the nth order included : `sh_first[0]` to `sh_first[n]`
/// 
/// ```
/// # use scilib::math::complex::Complex;
/// # use scilib::math::bessel::*;
/// 
/// let res = sh_first_array(Complex::from(13, 5), 3);
/// assert_eq!(res[3], Complex::from(5.127390110655217e-4, 2.5295761378174575e-4));
/// ```
pub fn sh_first_array<T>(z: T, n: usize) -> Vec<Complex> 
where T: Into<Complex> {
    let x: Complex = z.into();
    let sj_res = sj_array(x, n);
    let sy_res = sy_array(x, n);
    let mut sh_first = Vec::<Complex>::with_capacity(n + 1);
    for i in 0..sj_res.len() {
        sh_first.push(sj_res[i] + Complex::i() * sy_res[i]);
    }
    sh_first
}

/// # Second spherical Hankel function: h2
/// 
/// Compute the second kind of spherical hankel function.
/// * `z` - where the function is evaluated
/// * `n` - order evaluated
/// 
/// ```
/// # use scilib::math::complex::Complex;
/// # use scilib::math::bessel::*;
/// 
/// let res = sh_second(Complex::from(13, 5), 3);
/// assert_eq!(res, Complex::from(3.22144998903028, -8.645011398687984));
/// ```
pub fn sh_second<T>(z: T, n: usize) -> Complex 
where T: Into<Complex> {
    let x: Complex = z.into();
    sj(x, n) - Complex::i() * sy(x, n)
}

/// # Second spherical Hankel function (array): h2
/// 
/// Compute the second kind of spherical hankel function by reccurency to get an array of function.
/// * `z` - where the function is evaluated
/// * `n` - maximum order evaluated, return a vector from 0 to the nth order included : `sh_second[0]` to `sh_second[n]`
/// 
/// ```
/// # use scilib::math::complex::Complex;
/// # use scilib::math::bessel::*;
/// 
/// let res = sh_second_array(Complex::from(13, 5), 3);
/// assert_eq!(res[3], Complex::from(3.2214420145498694, -8.644990000503286));
/// ```
pub fn sh_second_array<T>(z: T, n: usize) -> Vec<Complex> 
where T: Into<Complex> {
    let x: Complex = z.into();
    let sj_res = sj_array(x, n);
    let sy_res = sy_array(x, n);
    let mut sh_second = Vec::<Complex>::with_capacity(n + 1);
    for i in 0..sj_res.len() {
        sh_second.push(sj_res[i] - Complex::i() * sy_res[i]);
    }
    sh_second
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
