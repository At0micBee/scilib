//!
//! # Distributions
//! 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

use std::f64::consts::{
    FRAC_1_PI,              // 1 / PI
    TAU                     // Tau constant
};

use super::{
    basic,                  // Using the basic functions
    series                  // Using the series tools
};          

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # Rayleigh distribution
/// 
/// ## Definition
/// The [Rayleigh distribution](https://en.wikipedia.org/wiki/Rayleigh_distribution) is a continuous
/// probability distribution for any $x>=0$. It is defined as:
/// $$
/// f_\sigma(x) = \frac{x}{\sigma^2} \exp\left( -\frac{x^2}{2\sigma^2} \right)
/// $$
/// 
/// ## Inputs
/// - `sigma`: the scale parameter of the distribution ($\sigma$)
/// - `x`: the value at which to evaluate the function ($x$).
/// 
/// Returns the density probability at a given point `x` for the `sigma`
/// scale parameter.
/// 
/// ## Example
/// ```
/// # use scilib::math::distribution::rayleigh;
/// let res = rayleigh(0.5, 0.9);
/// assert!((res - 0.7124353167010128).abs() < 1e-15);
/// ```
pub fn rayleigh(sigma: f64, x: f64) -> f64 {
    assert!(x >= 0.0);  // Rayleigh needs to be positive
    x * ( -0.5 * (x / sigma).powi(2) ).exp() / sigma.powi(2)
}

/* pub fn rayleigh_vec(sigma: f64, x: &[f64]) -> Vec<f64> {
    assert!(series::min_slice(x) >= 0.0);           // Rayleigh needs to be positive
    let frac_1_sigma2: f64 = 1.0 * sigma.powi(2);   // Pre-computing the factor

    x.iter().map(|v| v * frac_1_sigma2 * (-v.powi(2) * frac_1_sigma2 / 2.0)).collect()
} */

pub fn cauchy(gamma: f64, x0: f64, x: f64) -> f64 {
    FRAC_1_PI * gamma / ((x - x0).powi(2) + gamma.powi(2))
}

/* pub fn cauchy_vec(gamma: f64, x0: f64, x: &[f64]) -> Vec<f64> {
    let pi_gamma: f64 = FRAC_1_PI * gamma;
    let gamma2: f64 = gamma.powi(2);
    x.iter().map(|v| pi_gamma / ((v - x0).powi(2) + gamma2)).collect()
} */

pub fn laplace(b: f64, mu: f64, x: f64) -> f64 {
    (-(x - mu).abs() / b).exp() / (2.0 * b)
}

/* pub fn laplace_vec(b: f64, mu: f64, x: &[f64]) -> Vec<f64> {
    let frac_1_b: f64 = 1.0 / b;
    x.iter().map(|v| frac_1_b * (-(v - mu).abs() * frac_1_b).exp() / 2.0).collect()
} */

pub fn logistic(mu: f64, s: f64, x: f64) -> f64 {
    basic::sech((x - mu) / (2.0 * s)).powi(2) / (4.0 * s)
}

/* pub fn logistic_vec(mu: f64, s: f64, x: &[f64]) -> Vec<f64> {
    let frac_1_s2: f64 = 1.0 / (2.0 * s);
    x.iter().map(|v| frac_1_s2 * basic::sech(frac_1_s2 * (v - mu)).powi(2) / 2.0 ).collect()
} */

/// # Normal function
/// 
/// ## Definition
/// The [normal](https://en.wikipedia.org/wiki/Normal_distribution) is defined as:
/// $$
/// g(x) = \frac{1}{\sigma\sqrt{2\pi}}\cdot\exp\left( -\frac{(x - \mu)^2}{2\sigma^2} \right)
/// $$
/// 
/// ## Inputs
/// - `mu`: the expected value ($\mu$)
/// - `sigma`: the variance ($\sigma$)
/// - `x`: the value to evaluate ($x$)
/// 
/// Returns the value of the normal function with parameters
/// 
/// ## Example
/// ```
/// # use scilib::math::distribution::normal;
/// let res: f64 = normal(1.0, 2.0, 3.0);
/// assert!((res - 0.120985362259).abs() < 1e-12);
/// ```
pub fn normal(mu: f64, sigma: f64, x: f64) -> f64 {
    (-0.5 * ((x - mu) / sigma).powi(2)).exp() / (sigma * TAU.sqrt())
}

/// # Gaussian function
/// 
/// ## Definition
/// The [gaussian function](https://en.wikipedia.org/wiki/Gaussian_function) is defined as:
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
/// # use scilib::math::distribution::gaussian;
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
/// # use scilib::math::distribution::sigmoid;
/// let res1: f64 = sigmoid(-1.0);
/// let res2: f64 = sigmoid(0.0);
/// let res_comp: f64 = sigmoid(1.0);
/// assert!((res1 - 0.26894142136999).abs() < 1.0e-12);
/// assert_eq!(res2, 0.5);
/// assert_eq!(res1, 1.0 - res_comp);
/// ```
pub fn sigmoid(x: f64) -> f64 {
    1.0 / (1.0 + (-x).exp())
}
