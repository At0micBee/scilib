//!
//! # Basic math functions
//! 
//! This module provides access to many useful function that are not provided by the base Rust.
//! 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

use std::f64::consts::{     // Using std lib constants
    //FRAC_PI_2,              // Pi / 2
    FRAC_2_SQRT_PI,         // 2 / sqrt(Pi)
    TAU                     // Tau constant
};

use super::{                // Using parts from the crate
    super::constant,        // Calling scilib constants
    complex::Complex,       // Using Complex numbers
    polynomial::Bernoulli   // Bernoulli polynomials
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// Precision used for convergence
const PRECISION: f64 = 1.0e-12;

/// Stieltjes gamma computation precision
const STIELTJES_M: usize = 1_000_000;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # Sinus cardinal
/// 
/// The value of the [cardinal sinus](https://fr.wikipedia.org/wiki/Sinus_cardinal) is defined as:
/// $$
/// \mathrm{sinc}(x) = \frac{\sin(x)}{x}
/// $$
/// By convention, when $x = 0,~\mathrm{sinc}(x) = 1$.
/// 
/// ---
///
/// `x`: the value at which to evaluate the function ($x$).
/// 
/// Returns the sinc value of `x`.
/// 
/// ```
/// # use scilib::math::basic::sinc;
/// let x: f64 = sinc(0.0);
/// let y: f64 = sinc(1.2);
/// 
/// // Comparing to tabulated values
/// assert_eq!(x, 1.0);
/// assert!((y - 0.776699238306) < 1.0e-12);
/// ```
pub fn sinc(x: f64) -> f64 {
    if x == 0.0 {
        1.0
    } else {
        x.sin() / x
    }
}

/// # Binomial theorem
/// 
/// The [binomial theorem](https://en.wikipedia.org/wiki/Binomial_theorem) is defined as:
/// $$
/// \binom{n}{k} = \frac{n!}{k!(n - k)!}
/// $$
/// 
/// ---
///
/// `n` is the number of options ($n$) and `k` is the selection ($k$).
/// 
/// Returns `k` among `n`.
/// 
/// ```
/// # use scilib::math::basic::binomial;
/// let res: usize = binomial(4_usize, 2_usize);
/// 
/// assert_eq!(res, 6);
/// ```
/// 
/// ---
/// 
/// The implementation used here relies on the recurrence relation:
/// $$
/// \binom{n}{k} = \binom{n-1}{k} + \binom{n-1}{k-1}
/// $$
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
/// The [factorial function](https://en.wikipedia.org/wiki/Factorial) is defined as:
/// $$
/// n! = \prod_{i=1}^{n}i
/// $$
/// 
/// ---
/// 
/// `n` is the integer at which to evaluate the factorial ($n$).
/// 
/// Returns `n!`, the product of positive integers less or equal to `n`.
/// 
/// ```
/// # use scilib::math::basic::factorial;
/// let res: usize = factorial(5_usize);
/// 
/// assert_eq!(res, 120);
/// ```
pub fn factorial<T>(n: T) -> usize
where T: Into<usize> {
    (1..=n.into()).fold(1, |res, val| res * val)
}

/// # Stieltjes Gamma function
/// 
/// `n` is the order of the Stieltjes function to use.
/// `a` is the value at which to compute the function.
/// 
/// Returns the value of Gamma_n(a). At the moment, the results are only valid for the first few
/// orders, as the computation is very expansive. To get the basic Stieltjes coefficient, set `a=1`.
/// 
/// ```
/// # use scilib::math::basic::stieltjes;
/// let res1 = stieltjes(0, 1.0.into());
/// 
/// assert!((res1.re - 0.577215664).abs() <= 1e-6);
/// ```
pub fn stieltjes(n: usize, a: Complex) -> Complex {

    let mut res: Complex = - (a + STIELTJES_M as f64).ln().powi(n as i32 + 1) / (n as f64 + 1.0);

    for k in 0..STIELTJES_M {
        res += (a + k as f64).ln().powi(n as i32) / (a + k as f64);
    }

    res
}

/// # Hurwitz Zeta function
/// 
/// WARNING: still under development, results cannot be guarantee.
pub fn zeta<T, U>(s: T, a: U) -> Complex
where T: Into<f64>, U: Into<Complex> {

    // Conversions
    let a_c: Complex = a.into();
    let s_f: f64 = s.into();

    // If a is negative and even, we use Bernoulli
    if s_f == 0.0 || (s_f.is_sign_negative() && s_f % 2.0 == 0.0) {
        let ber: Bernoulli = Bernoulli::new(-s_f as usize + 1);
        return -ber.compute_complex(a_c) / (-s_f + 1.0);
    }

    let mut res: Complex = Complex::from(1.0 / (s_f - 1.0), 0.0);

    let mut n: usize = 0;
    let mut sign: f64 = -1.0;
    let mut div: f64;
    let mut term: Complex;

    'convergence: loop {

        if n >= 15 {
            break 'convergence;
        }

        sign *= -1.0;
        div = factorial(n) as f64;
        term = stieltjes(n, a_c);

        res += sign * term * (s_f - 1.0).powi(n as i32) / div;

        n += 1;
    }

    res
}

/// # Polylogarithm
pub fn li(s: f64, z: Complex) -> Complex {

    let mut n: usize = 1;
    let mut res_z: Complex = z;
    let mut div: Complex = (1.0_f64).powf(s).into();

    let mut term: Complex = res_z / div;

    if term.modulus() <= 1.0e-8 {
        return term;
    }

    let mut res: Complex = Complex::new();

    'convergence: loop {
        res += term;

        if (term / res).modulus() <= 1.0e-8 {
            break 'convergence;
        }

        n += 1;
        res_z *= z;
        div = (n as f64).powf(s).into();
        term = res_z / div;
    }

    res
}

/// # Gamma function
/// 
/// The [gamma function](https://en.wikipedia.org/wiki/Gamma_function) is a generalization of the factorial, and is defined as:
/// $$
/// \Gamma(z) = \int_{0}^{\infty}x^{z-1}\exp(-x)dx
/// $$
/// 
/// This function provides result for any real number, and returns the same result for integer as a factorial:
/// $$
/// \Gamma(n) = (n-1)!
/// $$
/// 
/// ---
/// 
/// `x` is the value to evaluate ($x$).
/// 
/// ```
/// # use scilib::math::basic::gamma;
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
pub fn gamma<T>(value: T) -> f64
where T: Into<f64> {

    let x: f64 = value.into();

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

/// # Euler Beta function
/// The [beta function](https://en.wikipedia.org/wiki/Beta_function) is an integral similar to the gamma function,
/// and it is defined as:
/// $$
/// B(x,y) = B(y,x) = \int_{0}^{1}t^{x-1}(1-t)^{y-1}dt
/// $$
/// 
/// ---
/// 
/// `x` and `y` are the points at which to evaluate the function ($x$, $y$).
/// 
/// The computation of the result is based on results from the gamma function. It is also possible to define
/// it with infinite series (or products), which could have some advantages, I'll investigate this possibility later.
/// 
/// ```
/// # use scilib::math::basic::beta;
/// let res: f64 = beta(1, 1.1);
/// let comp1: f64 = beta(3, 2);
/// let comp2: f64 = beta(2, 3);
/// 
/// assert!((res - 0.909090).abs() < 1.0e-5);
/// assert_eq!(comp1, comp2);
/// ```
/// 
/// ---
/// 
/// The current implementation relies on the relation:
/// $$
/// B(x,y) = \frac{\Gamma(x)\Gamma(y)}{\Gamma(x+y)}
/// $$
/// Which is easier to manage, but could be slower and slightly less precise.
/// Future updates will improve this function.
pub fn beta<T, U>(x: T, y: U) -> f64
where T: Into<f64> + Copy, U: Into<f64> + Copy {

    let t1: f64 = gamma(x);
    let t2: f64 = gamma(y);
    let b: f64 = gamma(x.into() + y.into());
    
    t1 * t2 / b
}

/// # Sigmoid function
/// The [sigmoid function](https://en.wikipedia.org/wiki/Sigmoid_function) is defined as:
/// $$
/// \sigma(x) = \frac{1}{1 + \exp(-x)} = 1 - \sigma(-x)
/// $$
/// 
/// ---
/// 
/// `x` is the value at which to evaluate the function ($x$).
/// 
/// Returns the value of the sigmoid function.
/// 
/// ```
/// # use scilib::math::basic::sigmoid;
/// let res1: f64 = sigmoid(-1.0);
/// let res2: f64 = sigmoid(0.0);
/// let res_comp: f64 = sigmoid(1.0);
/// 
/// assert!((res1 - 0.26894142136999).abs() < 1.0e-12);
/// assert_eq!(res2, 0.5);
/// assert_eq!(res1, 1.0 - res_comp);
/// ```
pub fn sigmoid(x: f64) -> f64 {
    1.0 / (1.0 + (-x).exp())
}

/// # Gaussian function
/// The [gaussian function](https://en.wikipedia.org/wiki/Gaussian_function) is defined as:
/// $$
/// g(x) = a\cdot\exp\left(-\frac{(b - x)^2}{2c^2}\right)
/// $$
/// 
/// ---
/// 
/// `a` is the amplitude ($a$), `b` is the center ($b$), 
/// `c` is the standard deviation ($c$),`x` is the value to evaluate ($x$).
/// 
/// Computes the value of the gaussian function with parameters a, b, c at x.
/// 
/// ```
/// # use scilib::math::basic::gaussian;
/// let res1: f64 = gaussian(1.0, 2.0, 3.0, 0.0);
/// 
/// assert_eq!(res1, 0.8007374029168081);
/// ```
pub fn gaussian(a: f64, b: f64, c: f64, x: f64) -> f64 {
    a * (-(x - b).powi(2) / (2.0 * c.powi(2))).exp()
}

/// # Normalized gaussian function
/// The [normalized gaussian function](https://en.wikipedia.org/wiki/Gaussian_function) is defined as:
/// $$
/// g(x) = \frac{1}{\sigma\sqrt{2\pi}}\cdot\exp\left( -\frac{(x - \mu)^2}{2\sigma^2} \right)
/// $$
/// 
/// ---
/// 
/// `lambda` is the expected value ($\mu$), `sigma` is the variance ($\sigma$).
/// 
/// ```
/// # use scilib::math::basic::gaussian_normed;
/// let res: f64 = gaussian_normed(1.0, 2.0, 3.0);
/// 
/// assert!((res - 0.120985362259).abs() < 1e-12);
/// ```
pub fn gaussian_normed(mu: f64, sigma: f64, x: f64) -> f64 {
    (1.0 / (sigma * (TAU.sqrt()))) * (-(x - mu).powi(2) / (2.0 * sigma.powi(2))).exp()
}

/// # Error function
/// The [error function](https://en.wikipedia.org/wiki/Error_function) is defined as:
/// $$
/// \mathrm{erf}(z) = \frac{2}{\sqrt{\pi}}\int_{0}^{z}\exp(-t^2)dt
/// $$
/// 
/// ---
/// 
/// `val` is the point at which to evaluate the function ($z$).
/// 
/// We define the error function for complex number.
/// 
/// WARNING: the erf function will soon become f64 only, moving the erf function for complex as a
/// complex function directly.
/// 
/// ```
/// # use scilib::math::complex::Complex;
/// # use scilib::math::basic::erf;
/// let r = erf(2.1);
/// let c = erf(Complex::from(-0.1, 0.7));
/// 
/// assert!((r.re - 0.997021).abs() < 1.0e-5);
/// assert!((c.re - -0.18297754).abs() < 1.0e-5 && (c.im - 0.92747498).abs() < 1.0e-5);
/// ```
pub fn erf<T>(val: T) -> Complex
where T: Into<Complex> {

    let x: Complex = val.into();

    let mut n: f64 = 0.0;               // Index of iteration
    let mut d1: f64 = 1.0;              // First div
    let mut d2: f64;                    // Second div
    let mut sg: f64 = 1.0;              // Sign of the term
    
    let mut term: Complex = x;          // Term at each iter
    let mut res: Complex = 0.0.into();  // Result

    // If the term is too small we exit
    if term.modulus() < PRECISION {
        return res;
    }

    'convergence: loop {
        res += term;

        // We exit when convergence reaches the precision
        if (term / res).modulus() < PRECISION {
            break 'convergence;
        }

        n += 1.0;
        sg *= -1.0;
        d1 *= n;
        d2 = 2.0 * n + 1.0;
        term = sg * x.powf(d2) / (d1 * d2);
    }

    FRAC_2_SQRT_PI * res
}

/// # Complementary error function
/// The [complementary error function](https://en.wikipedia.org/wiki/Error_function) is defined as:
/// $$
/// \mathrm{erfc}(z) = 1 - \mathrm{erf}(z)
/// $$
/// 
/// ---
/// 
/// `val` is the point at which to evaluate the function ($z$).
/// 
/// This function simply returns the complement of the error function, that is `1 - erf(z)`.
/// 
/// ```
/// # use scilib::math::complex::Complex;
/// # use scilib::math::basic::erfc;
/// let c = Complex::from(1.25, 0.3);
/// let res = erfc(c);
/// 
/// assert!((res.re - 0.0505570).abs() < 1.0e-5 && (res.im - -0.0663174).abs() < 1.0e-5);
/// ```
pub fn erfc<T>(val: T) -> Complex
where T: Into<Complex> {
    Complex::unity() - erf(val)
}

/// # Imaginary error function
/// The [imaginary error function](https://en.wikipedia.org/wiki/Error_function) is defined as:
/// $$
/// \mathrm{erfi}(z) = -i\cdot\mathrm{erf}(iz)
/// $$
/// 
/// ---
/// 
/// `val` is the point at which to evaluate the function ($z$).
/// 
/// This function is defined as `-i * erf(i*z)`.
/// 
/// ```
/// # use scilib::math::complex::Complex;
/// # use scilib::math::basic::erfi;
/// let c = Complex::from(0.07, -1.1);
/// let res = erfi(c);
/// 
/// assert!((res.re - 0.02349883).abs() < 1.0e-5 && (res.im - -0.88201955).abs() < 1.0e-5);
/// ```
pub fn erfi<T>(val: T) -> Complex
where T: Into<Complex> {
    -Complex::i() * erf(Complex::i() * val)
}

/// # Builds Pascal's triangle line
/// We can obtain any given line of the [pascal triangle](https://en.wikipedia.org/wiki/Pascal%27s_triangle) by using
/// the formula defined here:
/// $$
/// e_n^k = \binom{n}{k}
/// $$
/// Where $e_n^k$ is the element of the $n^{th}$ line and $i^{th}$ position.
/// 
/// ---
/// 
/// `n` the index of the line in the triangle.
/// 
/// Returns the line of index `n`.
/// 
/// The first position of the triangle is number as line 0, the rest follows from there.
/// 
/// ```
/// # use scilib::math::basic::pascal_triangle;
/// let l5: Vec<usize> = pascal_triangle(5);
/// 
/// assert_eq!(l5, vec![1, 5, 10, 10, 5, 1]);
/// ```
pub fn pascal_triangle(n: usize) -> Vec<usize> {
    let mut res: Vec<usize> = Vec::with_capacity(n);

    for k in 0..=n {
        res.push(binomial(n, k));
    }

    res
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
