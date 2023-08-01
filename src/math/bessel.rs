//!
//! # Bessel functions
//! 
//! The [Bessel functions](https://en.wikipedia.org/wiki/Bessel_function) are the solution to the Bessel differential equation.
//! There are multiple variants of these solutions, and this sub-module provides functions for all of them.
//! 
//! They are used to solve Bessel's differential equation:
//! $$
//! x^2\frac{d^2y}{dx^2} + x\frac{dy}{dx} + (x^2 - \alpha^2)y = 0
//! $$
//! 
//! ## First kind: $J_n$
//! 
//! The $J_n(x)$ function is the simplest solution to Bessel's equation. The most generalized case is for
//! a real index `n` and a complex number `x`, which is possible with this crate. Additionally, two methods are
//! provided for this equation, one for a integer order (named `j`) and one for a real order (named `jf`). The integer
//! order one is limited but is faster than its counterpart. The real order function is slower, but can compute any
//! order. Note that `jf` will attempt to fall back on `j` when it finds an integer order.
//! 
//! ```rust
//! # use num_complex::Complex64;
//! # use scilib::math::bessel::{ j_n, j_nu };
//! let f: f64 = 1.1;
//! let c = Complex64::new(-0.75, 3.0);
//! let res_i = j_n(1, f);                  // Faster for integer order
//! let res_f = j_nu(1.5, c);               // Would also work with 1.0
//! ```
//! 
//! ## Second kind: $Y_n$
//! 
//! Similar to the first kind, the $Y_n(x)$ equation are solution of Bessel's equation with a singularity at the origin.
//! The $Y$ function is itself based on the $J$ function. The function is undefined for any integer order, in which
//! case the limit has to be taken.
//! 
//! ```rust
//! # use num_complex::Complex64;
//! # use scilib::math::bessel::y;
//! let c = Complex64::new(2.0, -1.2);
//! let res_f = y(c, 1.5);              // Not a problem
//! let res_i = y(c, 1);                // The function takes the limit in this case
//! ```
//! 
//! ## Modified first kind: $I_n$
//! 
//! Also known as the hyperbolic Bessel function of the first kind. Its definition is similar to $J$, but lacks the
//! alternating $(-1)^k$ term in the sum.
//! 
//! ```rust
//! # use num_complex::Complex64;
//! # use scilib::math::bessel::i;
//! let c = Complex64::new(0.2, 1.0);
//! let res = i(c, -1.2);
//! ```
//! 
//! ## Modified second kind: $K_n$
//! 
//! Also known as the hyperbolic Bessel function of the second kind. Its definition is similar to $Y$, but lacks the
//! $cos(n\pi)$, and is normalized by $\pi/2$.
//! 
//! ```rust
//! # use num_complex::Complex64;
//! # use scilib::math::bessel::k;
//! let c = Complex64::new(0.0, 7.0);
//! let res = k(c, 0);
//! ```
//! 
//! ## Hankel functions: $H_n^{(1)}$ and $H_n^{(2)}$
//! 
//! Hankel functions are two linearly independent solutions to Bessel's equation.
//! 
//! ```rust
//! # use num_complex::Complex64;
//! # use scilib::math::bessel::{ hankel_first, hankel_second };
//! let c = Complex64::new(-0.3, 1.52);
//! let res_1 = hankel_first(c, -2.3);
//! let res_2 = hankel_second(c, -2.3);
//! ```
//! 
//! # Spherical Bessel functions
//! 
//! The spherical Bessel functions are used to solve the Helmholtz equation, which is:
//! $$
//! x^2\frac{d^2y}{dx^2} + 2x\frac{dy}{dx} + (x^2 - n(n+1))y = 0
//! $$
//! 
//! There is two kinds of spherical Bessel functions, plus the corresponding Hankel functions
//! 
//! ## First kind: $j_n$
//! First solution to the equation, the result is based on the $J_n$ Bessel equation.
//! 
//! ```rust
//! # use num_complex::Complex64;
//! # use scilib::math::bessel::{ sj, sj_array };
//! let c = Complex64::new(1.44, 2.22);
//! let res = sj(c, 2);         // Computes the single term
//! let res_a = sj_array(c, 2); // Computes all term up to 2
//! ```
//! 
//! ## Second kind: $y_n$
//! Second solution to the equation, the result of this one is based on the $Y_n$ Bessel function.
//! 
//! ```rust
//! # use num_complex::Complex64;
//! # use scilib::math::bessel::{ sy, sy_array };
//! let c = Complex64::new(2.53, -0.33);
//! let res = sy(c, 2);         // Computes the single term
//! let res_a = sy_array(c, 5); // Computes all term up to 5
//! ```
//! 
//! ## Hankel first kind: $h^{(1)}$
//! Follows the same definition as $H^{(1)}$, but using the spherical version of the Bessel equations.
//! 
//! ```rust
//! # use num_complex::Complex64;
//! # use scilib::math::bessel::{ sh_first, sh_first_array };
//! let c = Complex64::new(0.2253, 4.25);
//! let res = sh_first(c, 2);           // Computes the single term
//! let res_a = sh_first_array(c, 4);   // Computes all term up to 5
//! ```
//! 
//! ## Hankel second kind: $h^{(2)}$
//! Follows the same definition as $H^{(2)}$, but using the spherical version of the Bessel equations.
//! 
//! ```rust
//! # use num_complex::Complex64;
//! # use scilib::math::bessel::{ sh_second, sh_second_array };
//! let c = Complex64::new(-5.2, -0.356);
//! let res = sh_second(c, 2);          // Computes the single term
//! let res_a = sh_second_array(c, 7);  // Computes all term up to 7
//! ```
//! 
//! # Riccati-Bessel functions
//! 
//! The Riccati-Bessel functions are modified Bessel functions, and are solutions to the equation:
//! $$
//! x^2\frac{d^2y}{dx^2} + (x^2 - n(n + 1))y = 0
//! $$
//! 
//! ## First kind: $S_n$
//! Simplest solution to the differential equation, we base the computation off of $J_n(x)$.
//! 
//! ```rust
//! # use num_complex::Complex64;
//! # use scilib::math::bessel::riccati_s;
//! let c = Complex64::new(0.3, 1.22);
//! let res = riccati_s(c, 2);
//! ```
//! 
//! ## Second kind: $C_n$
//! Similar to $S_n$, but is based on the $Y_n(x)$ Bessel function.
//! 
//! ```rust
//! # use num_complex::Complex64;
//! # use scilib::math::bessel::riccati_c;
//! let c = Complex64::new(-4.2, 2.13);
//! let res = riccati_c(c, 5);
//! ```
//! 
//! ## Third kind: $\xi_n$
//! This solution makes use of the first kind of Hankel function $H^{(1)}$.
//! 
//! ```rust
//! # use num_complex::Complex64;
//! # use scilib::math::bessel::riccati_xi;
//! let c = Complex64::new(2.5, -0.25);
//! let res = riccati_xi(c, 4);
//! ```
//! 
//! ## Fourth kind: $\zeta_n$
//! This last solution makes use of the second kind of Hankel function $H^{(2)}$.
//! 
//! ```rust
//! # use num_complex::Complex64;
//! # use scilib::math::bessel::riccati_zeta;
//! let c = Complex64::new(-1.1, 8.2);
//! let res = riccati_zeta(c, 3);
//! ```

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

use std::f64::consts::{     // Using std lib constants
    PI,                     // Pi
    FRAC_PI_2               // Pi / 2
};

use super::basic;           // Basic functions

use num_complex::Complex64; // Using complex numbers from the num crate

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # Precision limit for Bessel computation
const PRECISION_CONVERGENCE: f64 = 1.0e-8;

const MAX_ITER_BESSEL: i32 = 200;

/// # Limit when computing Bessel Y
const DISTANCE_Y_LIM: f64 = 1e-4;

/// #
const PACK: usize = 50;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # $J_n$ Bessel function (integer order)
/// 
/// ## Definition
/// The $J_n$ Bessel represent the first kind of Bessel function, and are solution to Bessel differential equation,
/// and is defined as follows:
/// $$
/// J_n(x) = \sum_{k=0}^{\infty}\frac{(-1)^k}{k!(k+n)!}\left( \frac{x}{2} \right)^{2k+n}
/// $$
/// For any integer order `n`. Since for any real input, an integer order will also yield a real output,
/// this function is defined only for `f64`. For a complex input, see the `j_nu` function.
/// 
/// The bessel function depend on an infinite sum of terms; which we can't have.
/// The criterion chosen here is check each new term impacts the results significantly enough.
/// The default value selected in the program is defined by `const PRECISION_CONVERGENCE: f64 = 1.0e-8;`.
/// There is additionally a max limit on iteration, to avoid potential infinite loops in critical
/// values of the functions.
///
/// ## Inputs
/// - `n`: the order of the function ($n$)
/// - `x`: the value to evaluate ($x$)
/// 
/// Returns the value of the $n^{th}$ order of the Bessel $J_n$ function at $x$.
/// 
/// ## Example
/// ```
/// # use scilib::math::bessel::j_n;
/// // Computing some example values
/// let res_00 = j_n(0, 0.0);
/// let res_01 = j_n(0, 1.0);
/// let res_10 = j_n(1, 0.0);
/// let res_11 = j_n(1, 1.0);
/// let res_20 = j_n(2, 0.0);
/// let res_21 = j_n(2, 1.0);
/// let res = j_n(7, 5.2);
/// 
/// // Comparing to tabulated data
/// assert_eq!(res_00, 1.0);
/// assert!((res_01 - 0.7651976865).abs() < 1.0e-8);
/// assert_eq!(res_10, 0.0);
/// assert!((res_11 - 0.44005058).abs() < 1.0e-8);
/// assert_eq!(res_20, 0.0);
/// assert!((res_21 - 0.11490348).abs() < 1.0e-8);
/// assert!((res - 0.06544728).abs() < 1.0e-8);
/// 
/// // The function also handles negative orders
/// let pos1 = j_n(3, 3.2);
/// let neg1 = j_n(-3, 3.2);
/// let pos2 = j_n(6, 2.45);
/// let neg2 = j_n(-6, 2.45);
/// 
/// assert!(pos1 == -neg1);
/// assert!(pos2 == neg2);
/// ```
pub fn j_n(n: i32, x: f64) -> f64 {

    let np: i32 = n.abs();                                      // Getting the positive value of n
    let x2: f64 = x / 2.0;                                      // Halving x
    let mut k: i32 = 0;                                         // Order counter
    let mut d1: f64 = 1.0;                                      // First div
    let mut d2: f64 = basic::factorial(np as usize) as f64;     // Second div
    let mut term: f64 = x2.powi(np) / d2;                       // The term at each step
    let mut sum: f64 = f64::default();                          // The result of the operation

    // If the first term is already too small we exit directly
    if term.abs() < PRECISION_CONVERGENCE {
        return sum;
    }

    // Computing the terms of the infinite series
    'convergence: while k < MAX_ITER_BESSEL {

        sum += term;                                            // Updating the sum

        // If the changed compared to the final value is small we break
        if (term / sum).abs() < PRECISION_CONVERGENCE {
            break 'convergence;
        }

        k += 1;                                                 // Incrementing value, and flipping sign
        d1 *= -k as f64;                                        // Next value in the n! term
        d2 *= (np + k) as f64;                                  // Next value in the (n+k)! term
        term = x2.powi(np + 2 * k) / (d1 * d2);                 // Updating the sum term
    }

    // If n < 0 and is odd, we need to multiply by -1
    if n.is_negative() && np.rem_euclid(2) == 1 {
        -1.0 * sum
    } else {
        sum
    }
}

/// # $J_nu$ Bessel function (real order)
/// 
/// ## Definition
/// $$
/// J_\alpha(x) = \sum_{m=0}^{\infty}\frac{(-1)^m}{m!\Gamma(m+\alpha+1)}\left( \frac{x}{2} \right)^{2p+\alpha}
/// $$
/// 
/// Similar to the other J Bessel method, but this one allows the use of a real (float) index, rather
/// than an integer. This method is more costly to use than the other, and thus isn't recommended for
/// integer orders. Since non-integer order can yield complex values, the input and output of the
/// function are defined as so.
///
/// ## Inputs
/// - `order`: order of the function ($\alpha$)
/// - `x`: the value to evaluate ($x$)
/// 
/// Returns the value of the real order of the Bessel $J$ function at $x$.
/// 
/// ## Example
/// ```
/// # use num_complex::Complex64;
/// # use scilib::math::bessel::{ j_n, j_nu };
/// // This method allows the computation of real index for J
/// let res_pos: Complex64 = j_nu(2.5, 1.0.into());
/// let res_neg: Complex64 = j_nu(-1.75, 2.4.into());
/// assert!((res_pos.re - 0.04949681).abs() < 1.0e-4);
/// assert!((res_neg.re - 0.11990699).abs() < 1.0e-4);
/// 
/// // We can also check that the results are coherent with integers
/// let res_i: f64 = j_n(2, 0.75);
/// let res_r: Complex64 = j_nu(2.0, 0.75.into());
/// assert!((res_i - res_r.re).abs() < 1.0e-4);
/// 
/// // Because results are Complex, negative numbers are allowed
/// let neg: Complex64 = j_nu(2.3, (-0.75).into());
/// let expected: Complex64 = Complex64::new(0.0219887007, 0.030264850);
/// 
/// assert!((neg.re - expected.re).abs() < 1.0e-5 && (neg.im - expected.im).abs() < 1.0e-5);
/// 
/// // As for j, we can also use Complex numbers
/// let c: Complex64 = Complex64::new(1.2, 0.5);
/// let res: Complex64 = j_nu(1.5, c);
/// 
/// assert!((res.re - 0.3124202913).abs() < 1.0e-5 && (res.im - 0.1578998151) < 1.0e-5);
/// ```
pub fn j_nu(nu: f64, z: Complex64) -> Complex64 {

    // Required book-keeping to avoid NaNs on correct values
    if nu.is_sign_negative() && nu.fract() == 0.0 {             // If neg and integer
        if nu.rem_euclid(2.0) == 1.0 {                          // If odd we need to * -1.0
            return -1.0 * j_nu(nu.abs(), z);
        } else {                                                // If even it's just the positive value
            return j_nu(nu.abs(), z);
        }
    }

    let z2: Complex64 = z / 2.0;                                // Halving z
    let mut k: f64 = 0.0;                                       // Order counter
    let mut d1: f64 = 1.0;                                      // First div
    let mut d2: f64 = basic::gamma(nu + 1.0);                   // Second div
    let mut term: Complex64 = z2.powf(nu) / d2;                 // The term at each step
    let mut sum: Complex64 = Complex64::default();              // The result of the operation
    let mut counter: i32 = 0;                                   // Iteration counter
    
    // If the first term is already too small we exit directly
    if term.norm() < PRECISION_CONVERGENCE {
        return sum;
    }

    // Computing the terms of the infinite series
    'convergence: while counter < MAX_ITER_BESSEL {

        counter += 1;
        sum += term;

        // If the changed compared to the final value is small we break
        if (term / sum).norm() < PRECISION_CONVERGENCE {
            break 'convergence;
        }

        k += 1.0;                                               // Incrementing value
        d1 *= -k;                                               // Next value in the n! term
        d2 *= nu + k;                                           // Next value in the gamma(n+k+1) term
        term = z2.powf(k.mul_add(2.0, nu)) / (d1 * d2);
    }

    sum
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # Y Bessel function, real index
/// 
/// ## Definition
/// The $Y$ Bessel represent the second kind of Bessel function, and are solution to Bessel differential equation,
/// in the case of a singularity at the origin. It is defined as:
/// $$
/// Y_\alpha(x) = \frac{J_\alpha(x)\cos(\alpha\pi) - J_\alpha(x)}{\sin(\alpha\pi)}
/// $$
/// 
/// Because the function is not continuous for integer values of `n`, we need to compute the limit around these points.
/// We set the limit distance with `DISTANCE_Y_LIM`, compute the limit above and below the desired point and take the average.
/// We achieve precision under `1.0e-5` for non-integer`n`, and integer `n` using this approach.
/// 
/// ## Inputs
/// - `x`: the value to evaluate ($x$)
/// - `order`: the order of the function ($\alpha$)
/// 
/// Returns the value of the real order of the Bessel $Y$ function at $x$.
/// 
/// ## Example
/// ```
/// # use num_complex::Complex64;
/// # use scilib::math::bessel::y;
/// let res_pos = y(1.5, Complex64::new(1.0, 0.0));
/// let res_neg = y(-1.5, Complex64::new(1.0, 0.0));
/// 
/// assert!((res_pos.re - -1.10249557).abs() < 1.0e-5);
/// assert!((res_neg.re - -0.24029783).abs() < 1.0e-5);
/// 
/// // Values for integer n also works
/// let res_int_p = y(1.0, Complex64::new(0.5, 0.0));
/// let res_int_n = y(-1.0, Complex64::new(0.5, 0.0));
/// 
/// assert!((res_int_p.re - -1.47147239).abs() < 1.0e-5);
/// assert!((res_int_n.re - 1.47147239).abs() < 1.0e-5);
/// 
/// // We can compute negative value with Y, the result is complex
/// let res_neg = y(3.1, Complex64::new(-1.2, 0.0));
/// assert!((res_neg.re - 3.90596471).abs() < 1.0e-5);
/// assert!((res_neg.im - -1.32157214).abs() < 1.0e-5);
/// ```
pub fn y(nu: f64, z: Complex64) -> Complex64 {

    // If n is whole, we have to take the limit, otherwise it's direct
    if nu.fract() == 0.0 {
        (y(nu + DISTANCE_Y_LIM, z) + y(nu - DISTANCE_Y_LIM, z)) / 2.0
    } else {
        ((nu * PI).cos() * j_nu(nu, z) - j_nu(-nu, z)) / (nu * PI).sin()
    }
}

pub fn i_n(nu: f64, z: Complex64) -> Complex64 {

    todo!();
}

/// # $I$ modified Bessel function (real order)
/// 
/// ## Definition
/// The $I$ modified First Bessel function represent another kind of solution to the Bessel differential equation.
/// 
/// We use a definition of I based on an infinite series (similar to $J$). This way, we ensure good precision in
/// the computation, by following the definition:
/// $$
/// I_\alpha = i^{-\alpha}J_\alpha(ix)=\sum_{m=0}^{\infty}\frac{1}{m!\Gamma(m+\alpha+1)}\left( \frac{x}{2} \right)^{2m+\alpha}
/// $$
/// 
/// ## Inputs
/// - `x` is the value to evaluate ($x$)
/// - `order` the order of the function ($\alpha$)
/// 
/// Returns the value of the real order of the Bessel $I$ function at $x$.
/// 
/// ## Example
/// ```
/// # use num_complex::Complex64;
/// # use scilib::math::bessel::{i_nu, j_nu};
/// let res = i_nu(0.0, Complex64::new(1.2, 0.0));
/// assert!((res.re - 1.3937255841).abs() < 1.0e-6 && res.im == 0.0);
/// 
/// let c = Complex64::new(-1.2, 0.5);
/// let r2 = i_nu(-1.6, c);
/// assert!((r2.re - 0.549831309685).abs() < 1.0e-6);
/// assert!((r2.im - -0.123202230359).abs() < 1.0e-6);
/// 
/// // We can check that the values are coherent
/// let val = Complex64::new(3.2, -1.1);
/// let resi = i_nu(1.2, val);
/// let conf = j_nu(1.2, Complex64::i() * val) * Complex64::i().powf(-1.2);
/// assert!((resi.re - conf.re).abs() < 1.0e-14);
/// assert!((resi.im - conf.im).abs() < 1.0e-14);
/// ```
pub fn i_nu(nu: f64, z: Complex64) -> Complex64 {

    // Required book-keeping to avoid NaNs on correct values
    if nu.is_sign_negative() && nu.fract() == 0.0 {
        return i_nu(nu.abs(), z);
    }

    let z2: Complex64 = z / 2.0;                                // Halving z
    let mut k: f64 = 0.0;                                       // Order counter
    let mut d1: f64 = 1.0;                                      // First div
    let mut d2: f64 = basic::gamma(nu + 1.0);                   // Second div
    let mut term: Complex64 = z2.powf(nu) / d2;                 // The term at each step
    let mut sum: Complex64 = Complex64::default();              // The result of the operation
    let mut counter: i32 = 0;                                   // Iteration counter
    
    // If the first term is already too small we exit directly
    if term.norm() < PRECISION_CONVERGENCE {
        return sum;
    }

    // Computing the terms of the infinite series
    'convergence: while counter < MAX_ITER_BESSEL {
        
        counter += 1;
        sum += term;

        // If the changed compared to the final value is small we break
        if (term / sum).norm() < PRECISION_CONVERGENCE {
            break 'convergence;
        }

        k += 1.0;                                               // Incrementing value
        d1 *= k;                                                // Next value in the n! term
        d2 *= nu + k;                                           // Next value in the gamma(n+k+1) term
        term = z2.powf(k.mul_add(2.0, nu)) / (d1 * d2);
    }

    sum
}

/// # $K$ modified Bessel function
/// 
/// ## Definition
/// The $K$ modified Second Bessel function represent another kind of solution to the Bessel differential equation.
/// 
/// The definition of $K$ is similar to $Y$, but is based on $I$ and not $J$:
/// $$
/// K_\alpha(x) = \frac{\pi}{2}\frac{I_{-\alpha}(x) - I_\alpha(x)}{\sin(\alpha\pi)}
/// $$
/// 
/// ## Inputs
/// - `x` is the value to evaluate ($x$)
/// - `order` the order of the function ($\alpha$)
/// 
/// Returns the value of the real order of the Bessel $K$ function at $x$.
/// 
/// ## Example
/// ```
/// # use num_complex::Complex64;
/// # use scilib::math::bessel::k;
/// let c1 = Complex64::new(2.0, -1.0);
/// let res = k(-3.5, c1);
/// assert!((res.re - -0.321136273642).abs() < 1.0e-5);
/// assert!((res.im - 0.767517851045).abs() < 1.0e-5);
/// 
/// // Similar to Y, we take the limit for integer orders
/// let c2 = Complex64::new(-1.1, 0.6);
/// let res_i = k(1.0, c2);
/// assert!((res_i.re - -1.615394011004).abs() < 1.0e-5);
/// assert!((res_i.im - -2.105684608842).abs() < 1.0e-5);
/// ```
pub fn k(nu: f64, z: Complex64) -> Complex64 {

    // If n is whole, we have to take the limit, otherwise it's direct
    if nu.fract() == 0.0 {
        (k(nu + DISTANCE_Y_LIM, z) + k(nu - DISTANCE_Y_LIM, z)) / 2.0
    } else {
        (FRAC_PI_2 / (nu * PI).sin()) * (i_nu(-nu, z) - i_nu(nu, z))
    }
}

/// # First Hankel function: $H^{(1)}$
/// 
/// ## Definition
/// Computes the first kind of Hankel function, which is defined as:
/// $$
/// H_\alpha^{(1)}(x) = J_\alpha(x) + iY_\alpha(x)
/// $$
/// 
/// ## Inputs
/// - `x` is the value to evaluate ($x$)
/// - `order` the order of the function ($\alpha$)
/// 
/// Returns the value of the real order of the first hankel function at $x$.
/// 
/// ## Example
/// ```
/// # use num_complex::Complex64;
/// # use scilib::math::bessel::h1_nu;
/// let c1 = Complex64::new(-1.1, 2.3);
/// let r1 = h1_nu(1.0, c1);
/// assert!((r1.re - -0.011202683694).abs() < 1.0e-5);
/// assert!((r1.im - 0.055194689304).abs() < 1.0e-5);
/// 
/// let c2 = Complex64::new(5.2, -3.0);
/// let r2 = h1_nu(-2.35, c2);
/// assert!((r2.re - -4.280947776983).abs() < 1.0e-5);
/// assert!((r2.im - 3.212350280920).abs() < 1.0e-5);
/// ```
pub fn h1_nu(nu: f64, z: Complex64) -> Complex64 {
    j_nu(nu, z) + Complex64::i() * y(nu, z)
}

/// # Second Hankel function: $H^{(2)}$
/// 
/// ## Definition
/// Computes the first kind of Hankel function, which is defined as:
/// $$
/// H_\alpha^{(2)}(x) = J_\alpha(x) - iY_\alpha(x)
/// $$
/// 
/// ## Inputs
/// - `x` is the value to evaluate ($x$)
/// - `order` the order of the function ($\alpha$)
/// 
/// Returns the value of the real order of the second hankel function at $x$.
/// 
/// ## Example
/// ```
/// # use num_complex::Complex64;
/// # use scilib::math::bessel::h2_nu;
/// let c1 = Complex64::new(-1.1, 2.3);
/// let r1 = h2_nu(1.0, c1);
/// assert!((r1.re - -3.544212072638).abs() < 1.0e-5);
/// assert!((r1.im - 2.298353924176).abs() < 1.0e-5);
/// 
/// let c2 = Complex64::new(5.2, -3.0);
/// let r2 = h2_nu(-2.35, c2);
/// assert!((r2.re - -0.006818448986).abs() < 1.0e-5);
/// assert!((r2.im - -0.019369789719).abs() < 1.0e-5);
/// ```
pub fn h2_nu(nu: f64, z: Complex64) -> Complex64 {
    j_nu(nu, z) - Complex64::i() * y(nu, z)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # First spherical Bessel function: $j_n$ (real order)
/// 
/// ## Definition
/// We follow the definition based on the $J$ Bessel functions:
/// $$
/// j_n(x) = \sqrt{\frac{\pi}{2x}}J_{n+\frac{1}{2}}(x)
/// $$
/// 
/// ## Inputs
/// - `z`: where the function is evaluated ($x$)
/// - `n`: order evaluated ($n$)
/// 
/// Returns the first kind of spherical bessel function.
/// 
/// ## Example
/// ```
/// # use num_complex::Complex64;
/// # use scilib::math::bessel::sj_nu;
/// let res_f = sj_nu(2.25, Complex64::new(1.12, 0.0));
/// assert!((res_f.re - 0.049958287243).abs() < 1e-6);
/// 
/// // also works for integer order, and complex input
/// let res = sj_nu(3.0, Complex64::new(13.0, 5.0));
/// assert!((res.re - 1.6109825049200244).abs() < 1e-8);
/// assert!((res.im - -4.322382277910797).abs() < 1e-8);
/// ```
pub fn sj_nu(nu: f64, z: Complex64) -> Complex64 {
    (FRAC_PI_2 / z).sqrt() * j_nu(nu + 0.5, z)
}

fn sj_upward_recurrence<T>(z: T, n: usize) -> Vec<Complex64> 
where T: Into<Complex64> {

    let count = n + 1;
    let x: Complex64 = z.into();
    let mut jn = vec![Complex64::default(); count];

    jn[0] = x.sin() / x;
    jn[1] = x.sin() / x.powi(2) - x.cos() / x;

    for i in 1..=count - 2 {
        jn[i + 1] = (2 * i + 1) as f64 / x * jn[i] - jn[i - 1];
    }

    jn
}

fn sj_downward_recurrence<T>(z: T, nl: usize, nu: usize) -> Vec<Complex64> 
where T: Into<Complex64> {

    let count = nu - nl + 1;
    let x: Complex64 = z.into();
    let mut jn= vec![Complex64::default(); count + 2];

    jn[count + 1] = Complex64::default();
    jn[count] = Complex64::new(1.0, 0.0);

    for i in (1..=count).rev() {
        jn[i - 1] = (2 * (nl + i) + 1) as f64 / x * jn[i] - jn[i + 1];
    }
    jn.resize(count, Complex64::default());

    jn
}

/// # First spherical Bessel function (array): j
/// 
/// ## Inputs
/// Computes the first kind of spherical bessel function by recurrence to get an array of function.
/// - `z`: where the function is evaluated
/// - `n`: maximum order evaluated
/// 
/// Returns a vector from 0 to the nth order included : `sj[0]` to `sj[n]`.
/// 
/// ## Example
/// ```
/// # use num_complex::Complex64;
/// # use scilib::math::bessel::*;
/// let res1 = sj_array(Complex64::new(13.0, 5.0), 3);
/// assert_eq!(res1[0], Complex64::new(3.8248700377925635, 3.708547263317134));
/// assert_eq!(res1[1], Complex64::new(-3.357143112679857, 3.9747696875545517));
/// assert_eq!(res1[2], Complex64::new(-4.1924320794482135, -2.6499227040139104));
/// 
/// let res2 = sj_array(0.2, 25);
/// assert_eq!(res2[13], Complex64::new(3.835110596379198e-24, 0.0));
/// assert_eq!(res2[17], Complex64::new(5.910455642760406e-33, 0.0));
/// assert_eq!(res2[25], Complex64::new(1.125476749298975e-51, 0.0));
/// 
/// // We can also confirm the coherence of the results:
/// let val = Complex64::new(3.2, -1.1);
/// let resj = sj(val, 3);
/// let resa = sj_array(val, 3);
/// assert!((resj.re - resa[3].re).abs() < 1.0e-6);
/// assert!((resj.im - resa[3].im).abs() < 1.0e-6);
/// ```
pub fn sj_array<T>(z: T, n: usize) -> Vec<Complex64>
where T: Into<Complex64> {

    let x: Complex64 = z.into();

    if x.norm() > n as f64 / 2.0 {
        return sj_upward_recurrence(x, n);
    }

    let num_big_loop = n / PACK;
    let mut jn_all = Vec::<Vec<Complex64>>::new();

    for i in 0..num_big_loop {
        jn_all.push(sj_downward_recurrence(x, i * PACK, (i + 1) * PACK));
    }

    let rest = n % PACK;
    if rest != 0 {
        jn_all.push(sj_downward_recurrence(x, n - rest, n));
    }

    let mut jn = Vec::<Complex64>::with_capacity(n);
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

/// # Second spherical Bessel function: $y$
/// 
/// ## Definition
/// We follow the definition based on the $Y$ Bessel functions:
/// $$
/// y_n(x) = \sqrt{\frac{\pi}{2x}}Y_{n+\frac{1}{2}}(x)
/// $$
/// 
/// ## Inputs
/// - `z`: where the function is evaluated
/// - `n`: order evaluated
/// 
/// Returns the second kind of spherical bessel function.
/// 
/// ## Example
/// ```
/// # use num_complex::Complex64;
/// # use scilib::math::bessel::sy_nu;
/// let res = sy_nu(3.0, Complex64::new(13.0, 5.0));
/// assert!((res.re - 4.322629120777188).abs() < 1e-8);
/// assert!((res.im - 1.6104674841102558) < 1e-8);
/// ```
pub fn sy_nu(n: f64, z: Complex64) -> Complex64 {
    (FRAC_PI_2 / z).sqrt() * y(n + 0.5, z)
}

/// # Second spherical Bessel function (array): y
/// 
/// ## Inputs
/// Computes the second kind of spherical bessel function by recurrence to get an array of function.
/// - `z`: where the function is evaluated
/// - `n`: maximum order evaluated
/// 
/// Returns a vector from 0 to the nth order included : `sy[0]` to `sy[n]`.
/// 
/// ## Example
/// ```
/// # use num_complex::Complex64;
/// # use scilib::math::bessel::*;
/// let res1 = sy_array(Complex64::new(13.0, 5.0), 3);
/// assert_eq!(res1[0], Complex64::new(-3.7090299518957797, 3.8248379131516654));
/// assert_eq!(res1[1], Complex64::new(-3.9748349852610523, -3.356650136356049));
/// assert_eq!(res1[2], Complex64::new(2.6504303824601, -4.1922958025278));
/// let res2 = sy_array(0.2, 25);
/// 
/// assert_eq!(res2[13], Complex64::new(-4.82921204481494e22, 0.0));
/// assert_eq!(res2[17], Complex64::new(-2.417182573235861e31, 0.0));
/// assert_eq!(res2[25], Complex64::new(-8.711173815326792e49, 0.0));
/// 
/// // We can also confirm the coherence of the results:
/// let val = Complex64::new(3.2, -1.1);
/// let resj = sy(val, 3);
/// let resa = sy_array(val, 3);
/// assert!((resj.re - resa[3].re).abs() < 1.0e-6);
/// assert!((resj.im - resa[3].im).abs() < 1.0e-6);
/// ```
pub fn sy_array<T>(z: T, n: usize) -> Vec<Complex64>
where T: Into<Complex64> {

    let x: Complex64 = z.into();
    let y0 = -x.cos() / x;

    if n == 0 {
        return vec![y0];
    }

    let y1 = -x.cos() / (x * x) - x.sin() / x;

    if n == 1 {
        return vec![y0, y1];
    }

    let mut yn: Vec<Complex64> = vec![Complex64::default(); n + 1];
    yn[0] = y0;
    yn[1] = y1;

    for i in 1..=n - 1 {
        yn[i + 1] = ((2 * i + 1) as f64) / x * yn[i] - yn[i - 1];
    }

    yn
}

/// # First spherical Hankel function: $h^{(1)}$
/// 
/// ## Definition
/// The first spherical hankel function is defined as:
/// $$
/// h_n^{(1)}(x) = j_n(x) + iy_n(x)
/// $$
/// 
/// ## Inputs
/// Compute the first kind of spherical hankel function.
/// - `z`: where the function is evaluated ($x$)
/// - `n`: order evaluated ($n$)
/// 
/// # Example
/// ```
/// # use num_complex::Complex64;
/// # use scilib::math::bessel::sh1_nu;
/// let res = sh1_nu(3.0, Complex64::new(2.5, 0.2));
/// assert!((res.re - -0.056229554655).abs() < 1e-5);
/// assert!((res.im - -0.750570918335).abs() < 1e-5);
/// ```
pub fn sh1_nu(nu: f64, z: Complex64) -> Complex64 {
    sj_nu(nu, z) + Complex64::i() * sy_nu(nu, z)
}

/// # First spherical Hankel function (array): $h^{(1)}$
/// 
/// ## Inputs
/// Compute the first kind of spherical hankel function by recurrence to get an array of function.
/// - `z`: where the function is evaluated
/// - `n`: maximum order evaluated
/// 
/// Returns a vector from 0 to the nth order included : `sh_first[0]` to `sh_first[n]`.
/// 
/// # Example
/// ```
/// # use num_complex::Complex64;
/// # use scilib::math::bessel::*;
/// let res = sh_first_array(Complex64::new(13.0, 5.0), 3);
/// assert_eq!(res[3], Complex64::new(5.127390110655217e-4, 2.5295761378174575e-4));
/// 
/// // We can also confirm the coherence of the results:
/// let val = Complex64::new(3.2, -1.1);
/// let resj = sh_first(val, 3);
/// let resa = sh_first_array(val, 3);
/// assert!((resj.re - resa[3].re).abs() <= 1.0e-5);
/// assert!((resj.im - resa[3].im).abs() <= 1.0e-6);
/// ```
pub fn sh_first_array<T>(z: T, n: usize) -> Vec<Complex64> 
where T: Into<Complex64> {

    let x: Complex64 = z.into();
    let sj_res = sj_array(x, n);
    let sy_res = sy_array(x, n);
    let mut sh_first = Vec::<Complex64>::with_capacity(n + 1);

    for i in 0..sj_res.len() {
        sh_first.push(sj_res[i] + Complex64::i() * sy_res[i]);
    }

    sh_first
}

/// # Second spherical Hankel function: $h^{(2)}$
/// 
/// ## Definition
/// The second spherical hankel function is defined as:
/// $$
/// h_n^{(2)}(x) = j_n(x) - iy_n(x)
/// $$
/// 
/// ## Inputs
/// - `z`: where the function is evaluated ($x$)
/// - `n`: order evaluated ($n$)
/// 
/// Returns the second kind of spherical hankel function.
/// 
/// ## Example
/// ```
/// # use num_complex::Complex64;
/// # use scilib::math::bessel::sh2_nu;
/// let res = sh2_nu(3.0, Complex64::new(13.0, 5.0));
/// assert!((res.re - 3.22144998903028).abs() < 1e-8);
/// assert!((res.im + 8.645011398687984).abs() < 1e-8);
/// ```
pub fn sh2_nu(nu: f64, z: Complex64) -> Complex64 {
    sj_nu(nu, z) - Complex64::i() * sy_nu(nu, z)
}

/// # Second spherical Hankel function (array): $h^{(2)}$
/// 
/// ## Inputs
/// Compute the second kind of spherical hankel function by recurrence to get an array of function.
/// - `z`: where the function is evaluated
/// - `n`: maximum order evaluated
/// 
/// Returns a vector from 0 to the nth order included : `sh_second[0]` to `sh_second[n]`.
/// 
/// ## Example
/// ```
/// # use num_complex::Complex64;
/// # use scilib::math::bessel::*;
/// let res = sh_second_array(Complex64::new(13.0, 5.0), 3);
/// assert_eq!(res[3], Complex64::new(3.2214420145498694, -8.644990000503286));
/// 
/// // We can also confirm the coherence of the results:
/// let val = Complex64::new(3.2, -1.1);
/// let resj = sh_second(val, 3);
/// let resa = sh_second_array(val, 3);
/// assert!((resj.re - resa[3].re).abs() < 1.0e-6);
/// assert!((resj.im - resa[3].im).abs() < 1.0e-6);
/// ```
pub fn sh_second_array<T>(z: T, n: usize) -> Vec<Complex64> 
where T: Into<Complex64> {

    let x: Complex64 = z.into();
    let sj_res = sj_array(x, n);
    let sy_res = sy_array(x, n);
    let mut sh_second = Vec::<Complex64>::with_capacity(n + 1);

    for i in 0..sj_res.len() {
        sh_second.push(sj_res[i] - Complex64::i() * sy_res[i]);
    }
    
    sh_second
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # Riccati-Bessel function $S_n$
/// ## Definition
/// The $S_n$ Riccati-Bessel function is defined as:
/// $$
/// S_n(z) = \sqrt{\frac{\pi z}{2}}J_{n+\frac{1}{2}}(z)
/// $$
/// Where $J_{n+\frac{1}{2}}(z)$ is the Bessel function of the first kind.
/// 
/// ## Inputs
/// - `z` is the value to evaluate ($x$)
/// - `n` the order of the function ($n$)
/// 
/// ## Example
/// ```
/// # use std::f64::consts::FRAC_PI_2;
/// # use num_complex::Complex64;
/// # use scilib::math::bessel::{ j_nu, riccati_s };
/// let val: Complex64 = Complex64::new(2.0, -1.1);
/// let res: Complex64 = riccati_s(3, val);
/// let reference: Complex64 = j_nu(3.5, val) * (FRAC_PI_2 * val).sqrt();
/// 
/// // We can compare the results to a known value
/// assert!((res.re - -0.0417973).abs() < 1.0e-5);
/// assert!((res.im - -0.218255).abs() < 1.0e-5);
/// assert_eq!(res, reference);
/// ```
pub fn riccati_s(n: usize, z: Complex64) -> Complex64 {
    (FRAC_PI_2 * z).sqrt() * j_nu(n as f64 + 0.5, z)
}

/// # Riccati-Bessel function $C_n$
/// ## Definition
/// The $C_n$ Riccati-Bessel function is defined as:
/// $$
/// C_n(z) = -\sqrt{\frac{\pi z}{2}}Y_{n+\frac{1}{2}}(z)
/// $$
/// Where $Y_{n+\frac{1}{2}}(z)$ is the Bessel function of the second kind.
/// 
/// ## Inputs
/// - `z` is the value to evaluate ($x$)
/// - `n` the order of the function ($n$)
/// 
/// ## Example
/// ```
/// # use std::f64::consts::FRAC_PI_2;
/// # use num_complex::Complex64;
/// # use scilib::math::bessel::{ y, riccati_c };
/// let val: Complex64 = Complex64::new(-2.0, 1.2);
/// let res: Complex64 = riccati_c(3, val);
/// let reference: Complex64 = y(3.5, val) * -(FRAC_PI_2 * val).sqrt();
/// 
/// // We can compare the results to a known value
/// assert!((res.re - -0.8647593488).abs() <= 1.0e-5);
/// assert!((res.im - -1.19070243).abs() <= 1.0e-5);
/// assert_eq!(res, reference);
/// ```
pub fn riccati_c(n: usize, z: Complex64) -> Complex64 {
    -(FRAC_PI_2 * z).sqrt() * y(n as f64 + 0.5, z)
}

/// # Riccati-Bessel function $\xi_n$
/// ## Definition
/// The $\xi_n$ Riccati-Bessel function is defined as:
/// $$
/// \xi_n(z) = \sqrt{\frac{\pi z}{2}}H_{n+\frac{1}{2}}^{(1)}(z)
/// $$
/// Where $H_{n+\frac{1}{2}}^{(1)}(z)$ is the Hankel function of the first kind.
/// 
/// ## Inputs
/// - `z` is the value to evaluate ($x$)
/// - `n` the order of the function ($n$)
/// 
/// ## Example
/// ```
/// # use std::f64::consts::FRAC_PI_2;
/// # use num_complex::Complex64;
/// # use scilib::math::bessel::{ h1_nu, riccati_xi };
/// let val: Complex64 = Complex64::new(1.25, 2.2);
/// let res: Complex64 = riccati_xi(2, val);
/// let reference: Complex64 = h1_nu(2.5, val) * (FRAC_PI_2 * val).sqrt();
/// 
/// // We can compare the results to a known value
/// assert!((res.re - -0.273294).abs() <= 1.0e-5);
/// assert!((res.im - -0.0245767).abs() <= 1.0e-5);
/// assert_eq!(res, reference);
/// ```
pub fn riccati_xi(n: isize, z: Complex64) -> Complex64 {
    (FRAC_PI_2 * z).sqrt() * h1_nu(n as f64 + 0.5, z)
}

/// # Riccati-Bessel function $\zeta_n$
/// ## Definition
/// The $\zeta_n$ Riccati-Bessel function is defined as:
/// $$
/// \zeta_n(z) = \sqrt{\frac{\pi z}{2}}H_{n+\frac{1}{2}}^{(2)}(z)
/// $$
/// Where $H_{n+\frac{1}{2}}^{(2)}(z)$ is the Hankel function of the first kind.
/// 
/// ## Inputs
/// - `z` is the value to evaluate ($x$)
/// - `n` the order of the function ($n$)
/// 
/// ## Example
/// ```
/// # use std::f64::consts::FRAC_PI_2;
/// # use num_complex::Complex64;
/// # use scilib::math::bessel::{ h2_nu, riccati_zeta };
/// let val: Complex64 = Complex64::new(1.25, 2.2);
/// let res: Complex64 = riccati_zeta(-2, val);
/// let reference: Complex64 = h2_nu(-1.5, val) * (FRAC_PI_2 * val).sqrt();
/// 
/// // We can compare the results to a known value
/// assert!((res.re - -6.17727).abs() <= 1.0e-5);
/// assert!((res.im - -0.195811).abs() <= 1.0e-5);
/// assert_eq!(res, reference);
/// ```
pub fn riccati_zeta(n: isize, z: Complex64) -> Complex64 {
    (FRAC_PI_2 * z).sqrt() * h2_nu(n as f64 + 0.5, z)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
