//!
//! # Distributions
//! 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

use std::f64::consts::{
    FRAC_1_PI,              // 1 / PI
    TAU                     // Tau constant
};

use super::basic;           // Using the basic functions
          
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # Rayleigh distribution
/// 
/// ## Definition
/// The [Rayleigh distribution](https://en.wikipedia.org/wiki/Rayleigh_distribution) is a continuous
/// probability distribution for any $x\ge0$. It is defined as:
/// $$
/// f_\sigma(x) = \frac{x}{\sigma^2} \exp\left( -\frac{x^2}{2\sigma^2} \right)
/// $$
/// 
/// ## Inputs
/// - `sigma`: the scale parameter of the distribution ($\sigma$)
/// - `x`: the value at which to evaluate the function, must be $\ge0$ ($x$).
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
/// 
/// ![](https://raw.githubusercontent.com/At0micBee/scilib/dev/imgs/distribution/rayleigh.png)
pub fn rayleigh(sigma: f64, x: f64) -> f64 {
    assert!(x >= 0.0);  // Rayleigh needs to be positive
    x * ( -0.5 * (x / sigma).powi(2) ).exp() / sigma.powi(2)
}

/// # Cumulative function of the Rayleigh distribution
/// 
/// ## Definition
/// The cumulative function for Rayleigh is defined by:
/// $$
/// F_\sigma(x) = 1 - \exp\left(-\frac{x^2}{2\sigma^2}\right)
/// $$
/// 
/// ## Inputs
/// - `sigma`: the scale parameter of the distribution ($\sigma$)
/// - `x`: the value at which to evaluate the function, must be $\ge0$ ($x$).
/// 
/// Returns the cumulative density probability at a given point `x` for the `sigma`
/// scale parameter.
/// 
/// ## Example
/// ```
/// # use scilib::math::distribution::rayleigh_cumulative;
/// let res = rayleigh_cumulative(0.5, 0.9);
/// assert!((res - 0.8021013009163853).abs() < 1e-15);
/// ```
/// 
/// ![](https://raw.githubusercontent.com/At0micBee/scilib/dev/imgs/distribution/rayleigh_cumulative.png)
pub fn rayleigh_cumulative(sigma: f64, x: f64) -> f64 {
    assert!(x >= 0.0);  // Rayleigh needs to be positive
    1.0 - (-0.5 * (x / sigma).powi(2)).exp()
}

/// # Cauchy distribution
/// 
/// ## Definition
/// The [Cauchy distribution](https://en.wikipedia.org/wiki/Cauchy_distribution) is a continuous
/// probability distribution. It is defined as:
/// $$
/// f_\gamma(x_0, x) = \frac{\gamma}{\pi} \frac{1}{(x-x_0)^2 + \gamma^2}
/// $$
/// 
/// ## Inputs
/// - `gamma`: the scale parameter, HWHM ($\gamma$).
/// - `x0`: the location of the peak ($x_0$).
/// - `x`: the point at which to evaluate the function.
/// 
/// Returns the density probability at a given point `x` for the `gamma`
/// scale parameter centered on `x0`.
/// 
/// ## Example
/// ```
/// use scilib::math::distribution::cauchy;
/// let res = cauchy(0.5, -2.0, 0.9);
/// assert!((res - 0.018378168948255814).abs() < 1e-15)
/// ```
/// 
/// ![](https://raw.githubusercontent.com/At0micBee/scilib/dev/imgs/distribution/cauchy.png)
pub fn cauchy(gamma: f64, x0: f64, x: f64) -> f64 {
    FRAC_1_PI * gamma / ((x - x0).powi(2) + gamma.powi(2))
}

/// # Cumulative function of the Cauchy distribution
/// 
/// ## Definition
/// The cumulative function for Rayleigh is defined by:
/// $$
/// F_\gamma(x_0, x) = \frac{1}{\pi}\arctan\left( \frac{x-x_0}{\gamma} \right) + \frac{1}{2}
/// $$
/// 
/// ## Inputs
/// - `gamma`: the scale parameter, HWHM ($\gamma$).
/// - `x0`: the location of the peak ($x_0$).
/// - `x`: the point at which to evaluate the function.
/// 
/// Returns the cumulative density probability at a given point `x` for the `gamma`
/// scale parameter centered on `x0`.
/// 
/// ## Example
/// ```
/// use scilib::math::distribution::cauchy_cumulative;
/// let res = cauchy_cumulative(0.5, -2.0, 0.9);
/// assert!((res - 0.9456532942677374).abs() < 1e-15)
/// ```
/// 
/// ![](https://raw.githubusercontent.com/At0micBee/scilib/dev/imgs/distribution/cauchy_cumulative.png)
pub fn cauchy_cumulative(gamma: f64, x0: f64, x: f64) -> f64 {
    ((x - x0) / gamma).atan() * FRAC_1_PI + 0.5
}

/// # Laplace distribution
/// 
/// ## Definition
/// 
/// ## Inputs
/// 
/// ## Example
pub fn laplace(b: f64, mu: f64, x: f64) -> f64 {
    (-(x - mu).abs() / b).exp() / (2.0 * b)
}

/// # Logistic distribution
/// 
/// ## Definition
/// 
/// ## Inputs
/// 
/// ## Example
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
