//!
//! # Basic math functions
//! 
//! This module provides access to many useful function that are not provided by the base Rust.
//! 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

use std::f64::consts::FRAC_2_SQRT_PI;   // 2 / sqrt(Pi)

use super::super::constant;             // Calling scilib constants

use num_complex::Complex64;             // Using complex numbers from the num crate

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// Precision used for convergence
const PRECISION: f64 = 1.0e-12;

/// Maximum iteration count for convergence
const MAX_ITER: usize = 200;

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
/// assert!((y - 0.776699238306021958) < f64::EPSILON);
/// 
/// // Works for values close to zero
/// let sub: f64 = 1.0e-308;
/// assert!(sub.is_subnormal());
/// assert_eq!(sinc(f64::MIN_POSITIVE), 1.0);
/// assert_eq!(sinc(sub), 1.0);
/// ```
pub fn sinc(x: f64) -> f64 {
    if x != 0.0 {
        x.sin() / x
    } else {
        1.0
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
/// assert!((res - 4.082988165073596568).abs() <= 1e-15);
/// 
/// // Works for values close to zero
/// let sub: f64 = 1.0e-308;
/// assert!(sub.is_subnormal());
/// assert_eq!(coth(f64::MIN_POSITIVE), f64::INFINITY);
/// assert_eq!(coth(sub), f64::INFINITY);
/// ```
pub fn coth(x: f64) -> f64 {
    assert!(x != 0.0, "`coth` is undefined for 0!");    // Undefined for 0
    let e: f64 = (2.0 * x).exp();
    (e + 1.0) / (e - 1.0)
}

/// # Hyperbolic secant
/// 
/// ## Definition
/// The hyperbolic secant is defined as:
/// $$
/// \mathrm{sech}(x) = \frac{1}{\cosh(x)}
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
/// assert!((res - 0.552286154278204747).abs() <= 1e-15);
/// assert_eq!(sech(0.0), 1.0);
/// ```
pub fn sech(x: f64) -> f64 {
    x.cosh().recip()
}

/// # Hyperbolic cosecant
/// 
/// ## Definition
/// The hyperbolic cosecant is defined as:
/// $$
/// \mathrm{csch}(x) = \frac{1}{\sinh(x)}
/// $$
/// 
/// ## Inputs
/// - `x`: the value at which to evaluate the function ($x$).
/// 
/// Returns the csch value of `x`.
/// 
/// ## Example
/// ```
/// # use scilib::math::basic::csch;
/// let x: f64 = -1.2;
/// let res = csch(x);
/// assert_eq!(res, 1.0 / x.sinh());
/// assert!((res - -0.662487977194315480).abs() <= 1e-15);
/// ```
pub fn csch(x: f64) -> f64 {
    assert!(x != 0.0, "`csch` is undefined for 0!");    // Undefined for 0
    x.sinh().recip()
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

/// # Sigmoid function
/// 
/// ## Definition
/// The [sigmoid function](https://en.wikipedia.org/wiki/Sigmoid_function) is defined as:
/// $$
/// \sigma(x) = \frac{1}{1 + \exp(-x)} = 1 - \sigma(-x)
/// $$
/// 
/// ## Inputs
/// - `x`: the value at which to evaluate the function ($x$).
/// 
/// Returns the value of the sigmoid function.
/// 
/// ## Example
/// ```
/// # use scilib::math::basic::sigmoid;
/// let res1: f64 = sigmoid(-1.0);
/// let res2: f64 = sigmoid(0.0);
/// let res_comp: f64 = sigmoid(1.0);
/// assert!((res1 - 0.26894142136999).abs() < 1.0e-12);
/// assert_eq!(res2, 0.5);
/// assert_eq!(res1, 1.0 - res_comp);
/// ```
pub fn sigmoid(x: f64) -> f64 {
    (1.0 + (-x).exp()).recip()
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
    (1..=n.into()).product()
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
/// assert!((res_1 - 1.166711905198160345).abs() < 1.0e-5);
/// assert!((res_2 - -3.591387).abs() < 1.0e-5);
/// ```
pub fn gamma(x: f64) -> f64 {

    let mut n: usize = 0;                                       // Order counter
    let mut acc: f64 = 1.0;                                     // Accumulation term
    let mut term: f64;                                          // Term of each iteration

    // Computing the terms of the infinite series
    'convergence: loop {

        n += 1;                                                 // Incrementing counter
        term = (x / n as f64).exp() / (1.0 + x / n as f64);     // Computing term
        acc *= term;                                            // Accumulating

        //If the changes become too small, we stop
        if ((term - 1.0) / acc).abs() < PRECISION {
            break 'convergence;
        }
    }

    acc * (-x * constant::EULER_MASCHERONI).exp() / x
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
/// let res: f64 = beta(1.0, 1.1);
/// let comp1: f64 = beta(3.0, 2.0);
/// let comp2: f64 = beta(2.0, 3.0);
/// assert!((res - 0.909090).abs() < 1.0e-5);
/// assert_eq!(comp1, comp2);
/// ```
pub fn beta(x: f64, y: f64) -> f64 {
    gamma(x) * gamma(y) / gamma(x + y)
}

/// # Error function
/// 
/// ## Definition
/// The [error function](https://en.wikipedia.org/wiki/Error_function) is defined as:
/// $$
/// \mathrm{erf}(z) = \frac{2}{\sqrt{\pi}}\int_{0}^{z}\exp(-t^2)dt
/// $$
/// 
/// ## Inputs
/// - `val`: the point at which to evaluate the function ($z$)
/// 
/// Returns the error function value at $z$.
/// 
/// ## Example
/// ```
/// # use scilib::math::basic::erf;
/// let zero = erf(0.0);
/// let v1 = erf(-0.5);
/// let v2 = erf(2.1);
/// let large = erf(f64::INFINITY);
/// assert_eq!(zero, 0.0);
/// assert!((v1 - -0.520499877813046537).abs() < 1e-12);
/// assert!((v2 - 0.997020533343667014).abs() < 1e-12);
/// assert_eq!(large, 1.0);
/// ```
pub fn erf(x: f64) -> f64 {

    if x.abs() > 4.86431 {                      // Empirical limit of the series
        return x.signum();                      // Returning -1.0 or 1.0, depending on sign
    }

    let squared: f64 = - x.powi(2);             // x^2
    let mut term: f64;                          // The term for accumulation in the sum
    let mut sum: f64 = x;                       // The sum
    let mut prod: f64;                          // The product
    let mut n: usize = 0;                       // Iteration counter

    'convergence: while n < MAX_ITER {

        n += 1;                                 // Incrementing counter
        prod = 1.0;                             // Resetting product
        term = x / (2 * n + 1) as f64;          // Product for current n

        for k in 1..=n {
            prod *= squared / k as f64;         // Computing the product
        }

        term *= prod;                           // Updating term
        sum += term;                            // Updating sum

        // If the term has reached the given precision, we break
        if (term / sum).abs() < PRECISION {
            break 'convergence;
        }
    }
    
    FRAC_2_SQRT_PI * sum                        // Multiplying with the constant
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
/// # use scilib::math::basic::erfc;
/// let one = erfc(0.0);
/// let v1 = erfc(-0.5);
/// let v2 = erfc(2.1);
/// let large = erfc(f64::INFINITY);
/// assert_eq!(one, 1.0);
/// assert!((v1 - 1.520499877813046537).abs() < 1e-12);
/// assert!((v2 - 0.002979466656332985).abs() < 1e-12);
/// assert_eq!(large, 0.0);
/// ```
pub fn erfc(x: f64) -> f64 {
    1.0 - erf(x)
}

pub fn erf_complex(x: Complex64) -> Complex64 {

    let squared: Complex64 = - x.powi(2);             // x^2
    let mut term: Complex64;                          // The term for accumulation in the sum
    let mut sum: Complex64 = x;                       // The sum
    let mut prod: Complex64;                          // The product
    let mut n = 0;                                     // Iteration counter

    'convergence: while n < MAX_ITER {

        n += 1;                                 // Incrementing counter
        prod = Complex64::new(1.0, 0.0);                             // Resetting product
        term = x / (2 * n + 1) as f64;          // Product for current n

        for k in 1..=n {
            prod *= squared / k as f64;         // Computing the product
        }

        term *= prod;                           // Updating term
        sum += term;                            // Updating sum

        // If the term has reached the given precision, we break
        if (term / sum).norm() < PRECISION {
            break 'convergence;
        }
    }
    
    FRAC_2_SQRT_PI * sum                        // Multiplying with the constant
}


pub fn erfc_complex<T>(val: Complex64) -> Complex64 {
    Complex64::new(1.0, 0.0) - erf_complex(val)
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
pub fn erfi(val: Complex64) -> Complex64 {
    -Complex64::i() * erf_complex(Complex64::i() * val)
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

    // Vec with the right amount of space
    let mut res: Vec<usize> = Vec::with_capacity(n + 1);

    // We compute the first half
    for k in 0..=(n / 2) {
        res.push(binomial(n, k));
    }

    // Then we mirror the values
    for k in (1 + n / 2)..=n {
        res.push(res[n - k]);
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
