//!
//! # Basic math functions
//! 
//! This module provides access to many useful function that are not provided by the base Rust.
//! 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

use std::f64::consts::FRAC_2_SQRT_PI;   // 2 / sqrt(Pi)

use super::{                            // Using parts from the crate
    super::constant,                    // Calling scilib constants
    polynomial::Poly                    // Bernoulli polynomials
};

use num_complex::Complex64;             // Using complex numbers from the num crate

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// Precision used for convergence
const PRECISION: f64 = 1.0e-12;

/// Stieltjes gamma computation precision
const STIELTJES_M: usize = 1_000_000;

/// Sum initialization for ln(gamma)
const GAMMA_INIT: f64 = 0.999999999999997092;

/// Offset for the base value in ln(gamma)
const GAMMA_BASE_OFFSET: f64 = 5.2421875;

/// Offset to apply to ln(gamma)
const GAMMA_LN_OFFSET: f64 = 2.5066282746310005;

/// Coefficients to compute the ln gamma function.
const GAMMA_COEFS: [f64; 14] = [
    57.1562356658629235, -59.5979603554754912,
    14.1360979747417471, -0.491913816097620199,
    0.339946499848118887e-4, 0.465236289270485756e-4,
    -0.983744753048795646e-4, 0.158088703224912494e-3,
    -0.210264441724104883e-3, 0.217439618115212643e-3,
    -0.164318106536763890e-3, 0.844182239838527433e-4,
    -0.261908384015814087e-4, 0.368991826595316234e-5
];

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # Sinus cardinal
/// 
/// ## Definition
/// The value of the [cardinal sinus](https://fr.wikipedia.org/wiki/Sinus_cardinal) is defined as:
/// $$
/// \mathrm{sinc}(x) = \frac{\sin(x)}{x}
/// $$
/// By convention, when $x = 0,~\mathrm{sinc}(x) = 1$.
/// 
/// ## Inputs
/// - `x`: the value at which to evaluate the function ($x$).
/// 
/// Returns the sinc value of `x`.
/// 
/// ## Example
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

/// # Hyperbolic cotangent
/// 
/// ## Definition
/// The hyperbolic cotangent is defined as:
/// $$
/// \coth(x) = \frac{\cosh(x)}{\sinh(x)} = \frac{\exp(2x) + 1}{\exp(2x) - 1}
/// $$ 
/// 
/// ## Inputs
/// - `x`: the value at which to evaluate the function ($x$).
/// 
/// Returns the coth value of `x`.
/// 
/// ## Example
/// ```
/// # use scilib::math::basic::coth;
/// let x = 0.25;
/// let res = coth(x);
/// assert!((res - 4.082988165073596).abs() <= 1e-15);
/// ```
pub fn coth(x: f64) -> f64 {
    assert!(x != 0.0);              // Undefined for 0
    let e: f64 = (2.0 * x).exp();
    (e + 1.0) / (e - 1.0)
}

/// # Hyperbolic secant
/// 
/// ## Definition
/// The hyperbolic secant is defined as:
/// $$
/// \sech(x) = \frac{1}{\cosh(x)}
/// $$
/// 
/// ## Inputs
/// - `x`: the value at which to evaluate the function ($x$).
/// 
/// Returns the sech value of `x`.
/// 
/// ## Example
/// ```
/// # use scilib::math::basic::sech;
/// let x: f64 = 1.2;
/// let res = sech(x);
/// assert_eq!(res, 1.0 / x.cosh());
/// assert!((res - 0.5522861542782047).abs() <= 1e-15)
/// ```
pub fn sech(x: f64) -> f64 {
    1.0 / x.cosh()
}

/// # Hyperbolic cosecant
/// 
/// ## Definition
/// The hyperbolic cosecant is defined as:
/// $$
/// \sech(x) = \frac{1}{\sinh(x)}
/// $$
/// 
/// ## Inputs
/// - `x`: the value at which to evaluate the function ($x$).
/// 
/// Returns the sech value of `x`.
/// 
/// ## Example
/// ```
/// # use scilib::math::basic::csch;
/// let x: f64 = 1.2;
/// let res = csch(x);
/// assert_eq!(res, 1.0 / x.sinh());
/// assert!((res - 0.6624879771943154).abs() <= 1e-15)
/// ```
pub fn csch(x: f64) -> f64 {
    assert!(x != 0.0);              // Undefined for 0
    1.0 / x.sinh()
}

/// # Gaussian function
/// 
/// ## Definition
/// The [gaussian function](https://en.wikipedia.org/wiki/Gaussian_function) is a function defined as:
/// $$
/// g(x) = a\cdot\exp\left(-\frac{(b - x)^2}{2c^2}\right)
/// $$
/// 
/// ## Inputs
/// - `a`: the amplitude ($a$)
/// - `b`: the center ($b$)
/// - `c`: the standard deviation ($c$)
/// - `x`: the value to evaluate ($x$).
/// 
/// Returns the value of the gaussian function with parameters $a$, $b$, $c$ at $x$.
/// 
/// ## Example
/// ```
/// # use scilib::math::basic::gaussian;
/// let res1: f64 = gaussian(1.0, 2.0, 3.0, 0.0);
/// assert_eq!(res1, 0.8007374029168081);
/// ```
pub fn gaussian(a: f64, b: f64, c: f64, x: f64) -> f64 {
    a * (-(x - b).powi(2) / (2.0 * c.powi(2))).exp()
}

/// # Binomial theorem
/// 
/// ## Definition
/// The [binomial theorem](https://en.wikipedia.org/wiki/Binomial_theorem) is defined as:
/// $$
/// \binom{n}{k} = \frac{n!}{k!(n - k)!}
/// $$
/// 
/// The implementation used here relies on the recurrence relation:
/// $$
/// \binom{n}{k} = \binom{n-1}{k} + \binom{n-1}{k-1}
/// $$
/// This pushes overflow back for a few more terms.
/// 
/// ## Inputs
/// - `n`: the number of options ($n$) and `k` is the selection ($k$).
/// 
/// Returns `k` among `n`.
/// 
/// ## Example
/// ```
/// # use scilib::math::basic::binomial;
/// let res: usize = binomial(4_usize, 2_usize);
/// assert_eq!(res, 6);
/// ```
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

/// # Generalized binomial function
/// 
/// ## Definition
/// The [generalized binomial theorem](https://en.wikipedia.org/wiki/Binomial_theorem) is defined as:
/// $$
/// \binom{n}{k} = \frac{n^\overline{k}}{k!}
/// $$
/// 
/// ## Inputs
/// - `n`: the number of options ($n$), can be real.
/// - `k`: is the selection ($k$).
/// 
/// Returns `k` among `n`.
/// 
/// ## Example
/// ```
/// # use scilib::math::basic::binomial_reduced;
/// let res: f64 = binomial_reduced(4.2, 2_usize);
/// assert!((res - 6.72).abs() < 1.0e-10);
/// ```
pub fn binomial_reduced(n: f64, r: usize) -> f64 {

    falling_factorial(n, r) / factorial(r) as f64
}

/// # Factorial function
/// 
/// ## Definition
/// The [factorial function](https://en.wikipedia.org/wiki/Factorial) is defined as:
/// $$
/// n! = \prod_{i=1}^{n}i
/// $$
/// 
/// ## Inputs
/// - `n`: the integer at which to evaluate the factorial ($n$).
/// 
/// Returns `n!`, the product of positive integers less or equal to `n`.
/// 
/// ## Example
/// ```
/// # use scilib::math::basic::factorial;
/// let res: usize = factorial(5_usize);
/// assert_eq!(res, 120);
/// ```
pub fn factorial<T>(n: T) -> usize
where T: Into<usize> {
    (1..=n.into()).fold(1, |res, val| res * val)
}

/// # Rising factorial
/// 
/// ## Definition
/// The [rising factorial](https://en.wikipedia.org/wiki/Falling_and_rising_factorials) is a polynomial, which can
/// be computed directly with:
/// $$
/// x^{\overline{n}} = \prod_{k=0}^{n-1}(x+k)
/// $$
/// 
/// ## Inputs
/// - `x`: the value to pass to to the function ($x$).
/// - `n`: the integer at which to evaluate the rising factorial factorial ($n$).
/// 
/// Returns the value of the rising factorial $x^{\overline{n}}$.
/// 
/// ## Example
/// ```
/// # use scilib::math::basic::rising_factorial;
/// let res: f64 = rising_factorial(3.2, 5_usize);
/// assert!((res - 3119.80032).abs() < 1e-12);
/// ```
pub fn rising_factorial<T, U>(x: T, k: U) -> f64
where T: Into<f64>, U: Into<usize> {
    let z: f64 = x.into();
    (0..k.into()).fold(1.0, |res, val| res * (z + val as f64))
}

/// # Kummer function
/// 
/// ## Definition
/// [Kummer function](https://en.wikipedia.org/wiki/Confluent_hypergeometric_function), is a generalized
/// hypergeometric series defined by:
/// $$
/// M(a,b,z) = \sum_{n=0}^{\infty}\frac{a^{(n)}z^n}{b^{(n)}n!}
/// $$
/// 
/// ## Inputs
/// - `a`: first parameter, will be computed in the rising factorial ($a$)
/// - `b`: second parameter, will be computed in the rising factorial ($b$)
/// - `x`: the value to evaluate ($x$)
/// 
/// Returns the value of $M(a,b,x)$.
/// 
/// ## Example
/// ```
/// # use scilib::math::basic::kummer_function;
/// let res = kummer_function(0.2, 1.3, 2.0);
/// let expected = 1.577568150906;
/// assert!((res - expected).abs() < 1.0e-8);
/// ```
pub fn kummer_function(a: f64, b: f64, x: f64) -> f64 {

    // Preparing the base values
    let mut n: usize = 0;                   // Iteration counter
    let mut f_n: f64 = 1.0;                 // Factorial at iteration
    let mut x_n: f64 = 1.0;                 // x^n
    let mut a_n: f64 = 1.0;                 // a at iteration n
    let mut b_n: f64 = 1.0;                 // b at iteration n

    let mut term: f64 = 1.0;                // The term at each iteration
    let mut res: f64 = 0.0;                 // The sum

    'convergence: loop {

        if (term / res).abs() < PRECISION {
            break 'convergence;
        }
        res += term;

        // Exit condition
        if n > 50 {                         // We have a n!, it's gonna go down quickly
            break 'convergence;
        }

        n += 1;
        f_n *= n as f64;
        x_n *= x;
        a_n *= a + (n - 1) as f64;
        b_n *= b + (n - 1) as f64;
        term = a_n * x_n / (b_n * f_n);
    }

    res
}

/// # Falling factorial
/// 
/// ## Definition
/// The [falling factorial](https://en.wikipedia.org/wiki/Falling_and_rising_factorials) is a polynomial, which can
/// be computed directly with:
/// $$
/// x^{\underline{n}} = \prod_{k=0}^{n-1}(x-k)
/// $$
/// 
/// ## Inputs
/// - `x`: the value to pass to to the function ($x$).
/// - `n`: the integer at which to evaluate the falling factorial factorial ($n$).
/// 
/// Returns the value of the falling factorial $x^{\underline{n}}$.
/// 
/// ## Example
/// ```
/// # use scilib::math::basic::falling_factorial;
/// let res: f64 = falling_factorial(3.2, 5_usize);
/// assert!((res - -1.35168).abs() < 1e-12);
/// ```
pub fn falling_factorial<T, U>(x: T, k: U) -> f64 
where T: Into<f64>, U: Into<usize> {
    let z: f64 = x.into();
    (0..k.into()).fold(1.0, |res, val| res * (z - val as f64))
}

/// # Stieltjes Gamma function
/// 
/// ## Definition
/// Based on the Stieltjes coefficients, the [Stieltjes gamma function](https://en.wikipedia.org/wiki/Stieltjes_constants)
/// computes associated values, based on the formula:
/// $$
/// \gamma_n(x) = \lim_{m \to \infty} \sum_{k=0}^m \left( \frac{\ln(k+x)^n}{k+x} - \frac{\ln(m+x)^{n+1}}{n+1} \right)
/// $$
/// 
/// ## Inputs
/// - `n`: the order of the Stieltjes function to use.
/// - `a`: the value at which to compute the function.
/// 
/// Returns the value of Gamma_n(a). At the moment, the results are only valid for the first few
/// orders, as the computation is very expansive. To get the basic Stieltjes coefficient, set `a=1`.
/// 
/// ## Example
/// ```
/// # use scilib::math::basic::stieltjes;
/// let res1 = stieltjes(0, 1.0.into());
/// assert!((res1.re - 0.577215664).abs() <= 1e-6);
/// ```
pub fn stieltjes(n: usize, a: Complex64) -> Complex64 {

    let mut res: Complex64 = - (a + STIELTJES_M as f64).ln().powi(n as i32 + 1) / (n as f64 + 1.0);

    for k in 0..STIELTJES_M {
        res += (a + k as f64).ln().powi(n as i32) / (a + k as f64);
    }

    res
}

/// # Hurwitz Zeta function
/// 
/// WARNING: still under development, results cannot be guarantee.
pub fn zeta<T, U>(s: T, a: U) -> Complex64
where T: Into<f64>, U: Into<Complex64> {

    // Conversions
    let a_c: Complex64 = a.into();
    let s_f: f64 = s.into();

    // If a is negative and even, we use Bernoulli
    if s_f == 0.0 || (s_f.is_sign_negative() && s_f % 2.0 == 0.0) {
        let ber: Poly = Poly::bernoulli(-s_f as usize + 1);
        return -ber.compute_complex(a_c) / (-s_f + 1.0);
    }

    let mut res: Complex64 = Complex64::new(1.0 / (s_f - 1.0), 0.0);

    let mut n: usize = 0;
    let mut sign: f64 = -1.0;
    let mut div: f64;
    let mut term: Complex64;

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
/// 
/// ## Definition
/// The [polylogarithm](https://en.wikipedia.org/wiki/Polylogarithm) is defined by the infinite series:
/// $$
/// Li_s(z) = \sum_{k=1}^\infty \frac{z^k}{k^s}
/// $$
/// Because of these properties, the polylogarithm can be computed only for $|z|<= 1$.
/// 
/// ## Inputs
/// - `s`: the power of the divisor ($s$)
/// - `z` is the computed value ($z$).
/// 
/// Returns the value of the polylogarithm $Li_s(z)$.
/// 
/// ## Example
/// ```
/// # use num_complex::Complex64;
/// # use scilib::math::basic::li;
/// let val: Complex64 = Complex64::new(0.52, -0.55);
/// let res = li(1.35, val);
/// assert!((res.re - 0.38167313).abs() <= 1.0e-8);
/// assert!((res.im - -0.794472094).abs() <= 1.0e-8);
/// ```
pub fn li(s: f64, z: Complex64) -> Complex64 {

    let mut n: usize = 1;
    let mut res_z: Complex64 = z;
    let mut div: Complex64 = (1.0_f64).powf(s).into();

    let mut term: Complex64 = res_z / div;

    if term.norm() <= 1.0e-8 {
        return term;
    }

    let mut res: Complex64 = Complex64::default();

    'convergence: loop {
        res += term;

        if (term / res).norm() <= 1.0e-8 {
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
/// ## Definition
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
/// With the current computation scheme, we limit the precision of the computation in exchange for speed.
/// Typical values are achieve within a `1.0e-5` margin of error. Changing the method to another one
/// might grant some more speed and lower the error on the results.
/// 
/// ## Inputs
/// - `x`: the value to evaluate ($x$).
/// 
/// Returns the value of the gamma function.
/// 
/// ## Example
/// ```
/// # use scilib::math::basic::gamma;
/// let res_1: f64 = gamma(2.3);
/// let res_2: f64 = gamma(-0.45);
/// assert!((res_1 - 1.16671190).abs() < 1.0e-5);
/// assert!((res_2 - -3.591387).abs() < 1.0e-5);
/// ```
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

/// # Ln(Gamma) function
/// 
/// # Definition
/// Instead of computing the whole Gamma function, we can simply compute the natural log value of the function.
/// This helps to limit overflow of the value.
/// 
/// Lanczos, C. “A Precision Approximation of the Gamma Function.”
/// Journal of the Society for Industrial and Applied Mathematics: Series B, Numerical Analysis, vol. 1, 1964, pp. 86–96. JSTOR,
/// [http://www.jstor.org/stable/2949767](http://www.jstor.org/stable/2949767). Accessed 11 Apr. 2023.
/// 
/// ## Inputs
/// - `x`: the value to evaluate ($x$).
/// 
/// Returns the value of ln(gamma) function.
/// 
/// ## Example
/// ```
/// # use scilib::math::basic::gamma;
/// # use scilib::math::basic::ln_gamma;
/// let res_1: f64 = ln_gamma(1.2).exp();
/// let res_2: f64 = gamma(1.2);
/// assert!((res_1 - 0.9181687423997606).abs() < 1.0e-15);
/// assert!((res_2 - res_1).abs() < 1.0e-5);
/// ```
pub fn ln_gamma(x: f64) -> f64 {

    assert!(x > 0.0);                           // Only valid for x > 0

    let mut inc_x: f64 = x;                     // The incremented x in the loop

    let mut base: f64 = x + GAMMA_BASE_OFFSET;  // Base value
    base = (x + 0.5) * base.ln() - base;

    let res: f64 = GAMMA_COEFS.iter().fold(GAMMA_INIT, |sum, c| {
        inc_x += 1.0;
        sum + c / inc_x
    });

    base + (GAMMA_LN_OFFSET * res / x).ln()  // Final computation
}

/// # Regularized Gamma function P
/// 
/// ## Definition
/// The [regularized gamma functions](https://en.wikipedia.org/wiki/Incomplete_gamma_function) are functions
/// related to the gamma function. The first one, $P$, is the cumulative distribution function and is defined as:
/// $$
/// P(a, x) = \frac{\gamma(a, x)}{\Gamma(a)}
/// $$
/// Where $\gamma(a,x)$ is the incomplete lower gamma function and $\Gamma(a)$ is the gamma function.
/// 
/// ## Inputs
/// - `a`: the power in the series ($a$)
/// - `x`: the value to evaluate ($x$)
/// 
/// Returns the value of the regularized gamma function P.
/// 
/// ### Example
pub fn reg_gamma_p(a: f64, x: f64) -> f64 {

    // Checking that the range is correct;
    assert!(x >= 0.0, "x cannot be less than 0 in incomplete_gamma_p!");
    assert!(a > 0.0, "a cannot be less than or equal to 0 in incomplete_gamma_p!");

    // Calling the appropriate method to compute P
    if x < a + 1.0 {
        series_for_gamma(a, x)      // Series
    } else {
        1.0 - cf_for_gamma(a, x)    // Continued fraction
    }
}

/// # Regularized Gamma function Q
/// 
/// ## Definition
/// The [regularized gamma functions](https://en.wikipedia.org/wiki/Incomplete_gamma_function) are functions
/// related to the gamma function. The first one, $P$, is the cumulative distribution function and is defined as:
/// $$
/// Q(a, x) = \frac{\Gamma(a, x)}{\Gamma(a)} = 1 - P(a,x)
/// $$
/// Where $\Gamma(a,x)$ is the incomplete upper gamma function and $\Gamma(a)$ is the gamma function.
/// 
/// ## Inputs
/// - `a`: the power in the series ($a$)
/// - `x`: the value to evaluate ($x$)
/// 
/// Returns the value of the regularized gamma function Q.
/// 
/// ### Example
pub fn reg_gamma_q(a: f64, x: f64) -> f64 {

    // Checking that the range is correct;
    assert!(x >= 0.0, "x cannot be less than 0 in incomplete_gamma_p!");
    assert!(a > 0.0, "a cannot be less than or equal to 0 in incomplete_gamma_p!");

    // Calling the appropriate method to compute Q
    if x < a + 1.0 {
        1.0 - series_for_gamma(a, x)    // Series
    } else {
        cf_for_gamma(a, x)              // Continued fraction
    }
}

/// # Series representation for the regularized gamma function
fn series_for_gamma(a: f64, x: f64) -> f64 {

    // For x == 0 the result is always zero
    if x == 0.0 {
        return 0.0;
    }

    // Initializing the variables
    let lng: f64 = gamma(a).ln();
    let mut ap: f64 = a;
    let mut sum: f64 = 1.0 / a;
    let mut term: f64 = sum;

    // We compute the series term up to 100
    for _ in 0..100 {
        ap += 1.0;
        term *= x / ap;
        sum += term;

        // When we reach convergence we return the value
        if term.abs() < sum.abs() * PRECISION {
            return sum * (-x + a * x.ln() - lng).exp();
        }
    }
    
    // We return 0.0 if we couldn't compute the value
    0.0
}

/// # Continued fraction representation for the regularized gamma function
fn cf_for_gamma(a: f64, x: f64) -> f64 {

    // Initializing the variables
    let lng: f64 = gamma(a).ln();
    let mut an: f64;
    let mut term: f64;
    let mut b: f64 = x + 1.0 - a;
    let mut c: f64 = 1.0 / 1.0e-30;
    let mut d: f64 = 1.0 / b;
    let mut h: f64 = d;

    for i in 1..=100 {
        an = (-i * (i - a as isize)) as f64;
        b += 2.0;
        
        // Computing and correcting d if necessary
        d = an * d + b;
        if d.abs() < 1.0e-30 {
            d = 1.0e-30;
        }

        // Computing and correcting c if necessary
        c = b + an / c;
        if c.abs() < 1.0e-30 {
            c = 1.0e-30;
        }

        term = d * c;
        h *= term;

        // If we reach convergence we return
        if (term - 1.0).abs() < PRECISION {
            return h * (-x + a * x.ln() - lng).exp()
        }
    }

    // We return 0.0 if we couldn't compute the value
    0.0
}

/// # Euler Beta function
/// 
/// ## Definition
/// The [beta function](https://en.wikipedia.org/wiki/Beta_function) is an integral similar to the gamma function,
/// and it is defined as:
/// $$
/// B(x,y) = B(y,x) = \int_{0}^{1}t^{x-1}(1-t)^{y-1}dt
/// $$
/// 
/// The current implementation relies on the relation:
/// $$
/// B(x,y) = \frac{\Gamma(x)\Gamma(y)}{\Gamma(x+y)}
/// $$
/// Which is easier to manage, but could be slower and slightly less precise.
/// Future updates will improve this function.
/// 
/// ## Inputs
/// - `x` and `y` are the points at which to evaluate the function ($x$, $y$).
/// 
/// 
/// ## Example
/// ```
/// # use scilib::math::basic::beta;
/// let res: f64 = beta(1, 1.1);
/// let comp1: f64 = beta(3, 2);
/// let comp2: f64 = beta(2, 3);
/// assert!((res - 0.909090).abs() < 1.0e-5);
/// assert_eq!(comp1, comp2);
/// ```
pub fn beta<T, U>(x: T, y: U) -> f64
where T: Into<f64> + Copy, U: Into<f64> + Copy {

    let t1: f64 = gamma(x);
    let t2: f64 = gamma(y);
    let b: f64 = gamma(x.into() + y.into());
    
    t1 * t2 / b
}

/// # Error function
/// 
/// ## Definition
/// The [error function](https://en.wikipedia.org/wiki/Error_function) is defined as:
/// $$
/// \mathrm{erf}(z) = \frac{2}{\sqrt{\pi}}\int_{0}^{z}\exp(-t^2)dt
/// $$
/// 
/// We define the error function for complex number.
/// 
/// WARNING: the erf function will soon become f64 only, moving the erf function for complex as a
/// complex function directly.
/// 
/// ## Inputs
/// - `val`: the point at which to evaluate the function ($z$)
/// 
/// Returns the error function value at $z$.
/// 
/// ## Example
/// ```
/// # use num_complex::Complex64;
/// # use scilib::math::basic::erf;
/// let r = erf(2.1);
/// let c = erf(Complex64::new(-0.1, 0.7));
/// assert!((r.re - 0.997021).abs() < 1.0e-5);
/// assert!((c.re - -0.18297754).abs() < 1.0e-5 && (c.im - 0.92747498).abs() < 1.0e-5);
/// ```
pub fn erf<T>(val: T) -> Complex64
where T: Into<Complex64> {

    let x: Complex64 = val.into();

    let mut n: f64 = 0.0;                   // Index of iteration
    let mut d1: f64 = 1.0;                  // First div
    let mut d2: f64;                        // Second div
    let mut sg: f64 = 1.0;                  // Sign of the term
    
    let mut term: Complex64 = x;            // Term at each iter
    let mut res: Complex64 = 0.0.into();    // Result

    // If the term is too small we exit
    if term.norm() < PRECISION {
        return res;
    }

    'convergence: loop {
        res += term;

        // We exit when convergence reaches the precision
        if (term / res).norm() < PRECISION {
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
/// 
/// ## Definition
/// The [complementary error function](https://en.wikipedia.org/wiki/Error_function) is defined as:
/// $$
/// \mathrm{erfc}(z) = 1 - \mathrm{erf}(z)
/// $$
/// 
/// ## Inputs
/// - `val`: the point at which to evaluate the function ($z$)
/// 
/// Returns the complement of the error function, `1 - erf(z)`.
/// 
/// ## Example
/// ```
/// # use num_complex::Complex64;
/// # use scilib::math::basic::erfc;
/// let c = Complex64::new(1.25, 0.3);
/// let res = erfc(c);
/// assert!((res.re - 0.0505570).abs() < 1.0e-5 && (res.im - -0.0663174).abs() < 1.0e-5);
/// ```
pub fn erfc<T>(val: T) -> Complex64
where T: Into<Complex64> {
    Complex64::new(1.0, 0.0) - erf(val)
}

/// # Imaginary error function
/// 
/// ## Definition
/// The [imaginary error function](https://en.wikipedia.org/wiki/Error_function) is defined as:
/// $$
/// \mathrm{erfi}(z) = -i\cdot\mathrm{erf}(iz)
/// $$
/// 
/// ## Inputs
/// - `val`: the point at which to evaluate the function ($z$)
/// 
/// Returns the imaginary error function, $-i\mathrm{erf}(iz)$.
/// 
/// ## Example
/// ```
/// # use num_complex::Complex64;
/// # use scilib::math::basic::erfi;
/// let c = Complex64::new(0.07, -1.1);
/// let res = erfi(c);
/// assert!((res.re - 0.02349883).abs() < 1.0e-5 && (res.im - -0.88201955).abs() < 1.0e-5);
/// ```
pub fn erfi<T>(val: T) -> Complex64
where T: Into<Complex64> {
    -Complex64::i() * erf(Complex64::i() * val.into())
}

/// # Exponential integral
/// 
/// ## Definition
/// The [exponential integral](https://en.wikipedia.org/wiki/Exponential_integral) is a special case
/// of the incomplete gamma function an is defined as:
/// $$
/// E_{n}(x) = x^{n-1}\Gamma(1-n, x) = \int_1^\infty \frac{\exp(-xt)}{t^n}dt~\mathrm{for}~x>0,~n=0,1,...
/// $$
/// 
/// We take advantage of the fact that
/// $$
/// E_0(x) = \frac{\exp(-x)}{x}~\mathrm{and}~E_n(0) = \frac{1}{n-1}~\mathrm{for}~n>1
/// $$
/// 
/// ## Inputs
/// - `val`: The value to evaluate ($x$)
/// - `order`: The order of the exponential function to evaluate ($n$)
/// 
/// Returns the value of the $n$th order exponential integral for $x$.
/// 
/// ## Example
/// ```
/// # use scilib::math::basic::exp_int;
/// let res1 = exp_int(5.5, 3);
/// let res2 = exp_int(2.1, 0);
/// let res3 = exp_int(0.0, 5);
/// let res4 = exp_int(0.5, 1);
/// 
/// assert!((res1 - 0.0004987707).abs() < 1.0e-8);
/// assert!((res2 - 0.0583125848).abs() < 1.0e-8);
/// assert_eq!(res3, 0.25);
/// assert!((res4 - 0.5597735947).abs() < 1.0e-6);
/// ```
pub fn exp_int(val: f64, order: usize) -> f64 {

    // Checking the validity of the inputs
    assert!(
        !(val < 0.0 || (val == 0.0 && order <= 1)),
        "Invalid arguments in exponential integral!"
    );

    // We go through all the cases
    // If the order is zero it's easy
    if order == 0 {

        (-val).exp() / val
    
    // If the value is zero it's easy
    } else if val == 0.0 {

        1.0 / (order - 1) as f64
    
    // We compute the infinite fraction series
    } else if val > 1.0 {

        // Initializing variables
        let mut a: f64;
        let mut i_f64: f64;
        let mut term: f64;
        let mut b: f64 = val + order as f64;
        let mut c: f64 = 1.0 / 1.0e-30;
        let mut d: f64 = 1.0 / b;
        let mut h: f64 = d;
        let orderm1: f64 = (order - 1) as f64;

        // Iterating for the convergence
        for i in 1..=100 {
            i_f64 = i as f64;                   // Casting once to avoid doing it twice
            a = -i_f64 * (orderm1 + i_f64);
            b += 2.0;
            d = 1.0 / (a * d + b);
            c = b + a / c;
            term = c * d;
            h *= term;

            // If we have reached convergence we break
            if (term - 1.0).abs() < PRECISION {
                return h * (-val).exp();
            }

        }

        // If it doesn't, we return the last value for now
        // TODO: this should be changed to 0.0 probably
        h * (-val).exp()
    
    // For any other case we brute force the series
    } else {

        let orderm1: usize = order - 1;
        let orderm1_f64: f64 = orderm1 as f64;

        // Initializing the first term
        let mut res: f64 = match orderm1 {
            0 => - val.ln() - constant::EULER_MASCHERONI,
            _ => 1.0 / orderm1_f64
        };

        // Initializing the rest of the variables
        let mut i_f64: f64;
        let mut term: f64;
        let mut psi: f64;
        let mut f: f64 = 1.0;

        // Iterating for the convergence
        for i in 1..=100 {
            i_f64 = i as f64;
            f *= - val / i_f64;

            if i != orderm1 {
                term = -f / (i_f64 - orderm1_f64);
            } else {
                psi = -constant::EULER_MASCHERONI;

                for j in 1..orderm1 {
                    psi += 1.0 / j as f64;
                }

                term = f * (- val.ln() + psi);
            }

            // Updating the result
            res += term;

            // If we have reached our precision we return
            if term.abs() < res.abs() * PRECISION {
                return res;
            }
        }

        // If it doesn't, we return the last value for now
        // TODO: this should be changed to 0.0 probably
        res
    }
}

/// # Principal exponential integral
/// 
/// ## Definition
/// The [principal exponential integral](https://mathworld.wolfram.com/ExponentialIntegral.html) is defined as:
/// $$
/// \mathrm{Ei}(x) = \int_{-\infty}^{x}\frac{\exp(t)}{t}dt~\mathrm{for}~x>0
/// $$
/// 
/// ## Inputs
/// - `val`: the value to evaluate ($x$)
/// 
/// Returns the value of the principal exponential integral for $x$.
/// 
/// ## Example
/// ```
/// # use scilib::math::basic::exp_int_i;
/// let res1 = exp_int_i(3.2);
/// let res2 = exp_int_i(1.0e-9);
/// 
/// assert!((res1 - 11.3673026569).abs() < 1.0e-8);
/// assert!((res2 - -20.1460501710).abs() < 1.0e-8);
/// ```
pub fn exp_int_i(val: f64) -> f64 {

    assert!(val >= 0.0, "Invalid arguments in principal exponential integral!");

    // We check for special cases
    if val < 1.0e-30 {

        val.ln() + constant::EULER_MASCHERONI

    // If we can we use the power series
    } else if val <= - PRECISION.ln() {

        let mut i_f64: f64;
        let mut ct: f64;
        let mut sum: f64 = 0.0;
        let mut fact: f64 = 1.0;

        // We compute the convergence
        for i in 1..=100 {
            i_f64 = i as f64;
            fact *= val / i_f64;
            ct = fact / i_f64;
            sum += ct;

            // If convergence is reached
            if ct < PRECISION * sum {
                return sum + val.ln() + constant::EULER_MASCHERONI;
            }
        }

        // If it doesn't, we return the last value for now
        // TODO: this should be changed to 0.0 probably
        sum + val.ln() + constant::EULER_MASCHERONI

    // Else we brute force the series
    } else {

        // Initializing the variables
        let mut i_f64: f64;
        let mut ct: f64;
        let mut sum: f64 = 0.0;
        let mut fact: f64 = 1.0;

        // We compute the convergence
        for i in 1..=100 {
            i_f64 = i as f64;
            ct = fact;
            fact *= i_f64 / val;

            // We check the evolution of the control term
            // If we're under the required precision, we return
            if fact < PRECISION {
                break;      // TODO: this should be changed to 0.0 probably
        
            // If It keeps getting smaller we keep going
            } else if fact < ct {
                sum += ct;

            // If it stop diminishing we return
            } else {
                sum -= ct;
                break;      // TODO: this should be changed to 0.0 probably
            }
        }

        // We return the value
        val.exp() * (1.0 + sum) / val
    }
}

/// # Builds Pascal's triangle line
/// 
/// ## Definition
/// We can obtain any given line of the [pascal triangle](https://en.wikipedia.org/wiki/Pascal%27s_triangle) by using
/// the formula defined here:
/// $$
/// e_n^k = \binom{n}{k}
/// $$
/// Where $e_n^k$ is the element of the $n^{th}$ line and $i^{th}$ position.
/// 
/// The first position of the triangle is number as line 0, the rest follows from there.
/// 
/// ## Inputs
/// - `n`: the index of the line in the triangle.
/// 
/// Returns the line at index `n`.
/// 
/// ## Example
/// ```
/// # use scilib::math::basic::pascal_triangle;
/// let l5: Vec<usize> = pascal_triangle(5);
/// assert_eq!(l5, vec![1, 5, 10, 10, 5, 1]);
/// ```
pub fn pascal_triangle(n: usize) -> Vec<usize> {
    let mut res: Vec<usize> = Vec::with_capacity(n);

    for k in 0..=n {
        res.push(binomial(n, k));
    }

    res
}

/// # Levi-Civita symbol
/// 
/// ## Definition
/// $$
/// \epsilon_{ijk...l}
/// $$
/// Computes the result of the $n$ dimensional Levi-Civita symbol. Returns 1 if the indexes are an even permutation
/// and -1 if odd permutation. If there is any repetition it returns 0.
/// 
/// The function does not perform a check to ensure that all numbers are account for.
/// 
/// ## Inputs
/// - `val`: the list of the indexes to evaluate
/// 
/// Returns the value of the Levi-Civita symbol.
/// 
/// ## Example
/// ```
/// # use scilib::math::basic::levi_civita;
/// let pos = levi_civita(vec![2, 3, 4, 5, 1]);
/// let neg = levi_civita(vec![3, 2, 1]);
/// let z = levi_civita(vec![1, 1, 2, 3, 4, 5, 6]);
/// assert_eq!(pos, 1);
/// assert_eq!(neg, -1);
/// assert_eq!(z, 0);
/// ```
pub fn levi_civita(val: Vec<isize>) -> isize {

    let mut res: isize = 1;

    for idx1 in 0..val.len() {
        for idx2 in (idx1+1)..val.len() {
            res *= (val[idx2] - val[idx1]).signum();
        }
    }
    
    res
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
