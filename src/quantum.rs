//!
//! # Methods for quantum mechanics
//! 
//! This crate is dedicated to providing practical methods to solve quantum mechanics problems. In addition to
//! simple tools such as computing the $n,~l,~m$ numbers, we also provide the requirements to solve the 
//! wave-function:
//! $$
//! \Psi_{n,l,m}(r, \theta, \phi) = R_n^l(r)Y_l^m(\theta, \phi)
//! $$
//! 
//! Where $R_n^l(r)$ is the radial wave-function and $Y_l^m(\theta, \phi)$ is the angular wave-function (or spherical harmonics function).
//! Which have the respective forms:
//! $$
//! R_n^l(r) = \sqrt{\left( \frac{2}{na_B} \right)^3 \frac{(n-l-1)!}{2n\cdot(n+l)!}} \cdot \left( \frac{2r}{na_B} \right)^l
//! L_{n-l-1}^{2l+1}\left( \frac{2r}{na_B} \right) \exp\left( -\frac{r}{na_B} \right)
//! $$
//! $$
//! Y_l^m(\theta, \phi) = (-1)^m \sqrt{\frac{(2l+1)}{4\pi}\frac{(l-m)!}{(l+m)!}} P_l^m(\cos(\theta)) \exp(im\phi)
//! $$
//! With $L_n^m$ the Laguerre polynomials and $P_n^m$ the Legendre polynomials.

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

use std::f64::consts::{     // Using std lib constants
    PI                      // Pi
};

use crate::{                // Calling other modules
    math::{                 // Math crate
        basic,              // Basic functions
        complex::Complex,   // Using Complex numbers
        polynomial          // Special polynomials
    },
    constant as cst         // Calling scilib constants
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # Azimuthal quantum number $l$
/// The [azimuthal quantum number](https://en.wikipedia.org/wiki/Azimuthal_quantum_number) is
/// computed based on the value of the principal quantum number $n$.
/// It takes the integer values:
/// $$
/// 0 \le l \le n-1
/// $$
/// 
/// ```
/// # use scilib::quantum::get_l;
/// let l = get_l(3_usize);
/// assert_eq!(l, vec![0, 1, 2]);
/// ```
pub fn get_l(n: usize) -> Vec<usize> {
    (0..n).collect()
}

/// # Magnetic quantum number $m$
/// The [magnetic quantum number](https://en.wikipedia.org/wiki/Magnetic_quantum_number) is
/// computed based on the value of the azimuthal quantum number $l$.
/// It takes the integer values:
/// $$
/// -l \le m \le l
/// $$
/// 
/// ```
/// # use scilib::quantum::get_m;
/// let m = get_m(2);
/// assert_eq!(m, vec![-2, -1, 0, 1, 2]);
/// ```
pub fn get_m(l: i32) -> Vec<i32>{
    (-l..=l).collect()
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # Radial wave-function
/// Computes the result of the radial wave-function, as defined by:
/// $$
/// R_n^l(r) = \sqrt{\left( \frac{2}{na_B} \right)^3 \frac{(n-l-1)!}{2n\cdot(n+l)!}} \cdot \left( \frac{2r}{na_B} \right)^l
/// L_{n-l-1}^{2l+1}\left( \frac{2r}{na_B} \right) \exp\left( -\frac{r}{na_B} \right)
/// $$
/// With $L_n^m$ the Laguerre polynomials.
/// 
/// This function derives the normalization factor and the associated Laguerre polynomial
/// to compute any wave function. We provide the principal quantum number `n` ($n$) and the
/// azimuthal quantum number `l` ($l$), as well as the radius at which to compute the solution.
/// 
/// ```
/// # use scilib::quantum::radial_wavefunction;
/// // Computing the Rnl for n=2, l=1
/// let res = radial_wavefunction(2, 1, 1.3e-12);
/// ```
pub fn radial_wavefunction(n: usize, l: usize, r: f64) -> f64 {

    // Pre-computing
    let factor: f64 = r / (n as f64 * cst::A_0);

    // Computing the norm of the function
    let mut norm: f64 = (2.0 / (n as f64 * cst::A_0)).powi(3);
    norm *= basic::factorial(n - l - 1) as f64 / (2 * n * basic::factorial(n + l)) as f64;

    // Computing the term associated to the Laguerre polynomial
    let poly: f64 = polynomial::Laguerre::new(n - l - 1, 2 * l as i32 + 1).compute(2.0 * factor);

    // Finishing computation and returning
    (2.0 * factor).powi(l as i32) * norm.sqrt() * poly * (-factor).exp()
}

/// # Radial wave-function for vectors
/// 
/// Similar to `radial_wavefunction`, but computes the resulting values for a vector.
/// This has the advantages to speed things up when a lot of values are required,
/// as the norm doesn't have to be computed at each pass, as well as the polynomial.
/// 
/// ```
/// # use scilib::range;
/// # use scilib::quantum::radial_vec;
/// // Computing the Rnl for n=2, l=1
/// let x: Vec<f64> = range::linear(0, 1e-9, 500);
/// let res = radial_vec(2, 1, &x);
/// ```
pub fn radial_vec(n: usize, l: usize, r: &[f64]) -> Vec<f64> {

    // Preparing the poly
    let poly = polynomial::Laguerre::new(n - l - 1, 2 * l as i32 + 1);

    // Computing the norm of the function
    let mut norm: f64 = (2.0 / (n as f64 * cst::A_0)).powi(3);
    norm *= basic::factorial(n - l - 1) as f64 / (2 * n * basic::factorial(n + l)) as f64;
    norm = norm.sqrt();

    // Div part of the factor
    let div: f64 = 1.0 / (n as f64 * cst::A_0);

    // We initialize the values of the result with the polynomial
    let mut res: Vec<f64> = r.iter().map(|rad| poly.compute(2.0 * rad * div)).collect();

    // Computing the other parts of the function
    for (elem, rad) in res.iter_mut().zip(r) {
        *elem *= (2.0 * rad * div).powi(l as i32) * norm * (-rad * div).exp();
    }
    
    res
}

/// # Spherical harmonics
/// Provides the solution to the angular $Y_l^m(\theta, \phi)$ wave-function:
/// $$
/// Y_l^m(\theta, \phi) = (-1)^m \sqrt{\frac{(2l+1)}{4\pi}\frac{(l-m)!}{(l+m)!}} P_l^m(\cos(\theta)) \exp(im\phi)
/// $$
/// With $P_n^m$ the Legendre polynomials.
/// 
/// For any `l` ($l$) the azimuthal quantum number and `m` ($m$) the magnetic quantum number. The equation
/// produces the solution for a given set of angles $\theta$ and $\phi$.
/// 
/// ```
/// # use scilib::quantum::spherical_harmonics;
/// // Computing the Ylm for l=2, m=1 at theta = 0.2rad and phi = -0.7rad
/// let res = spherical_harmonics(2, 1, 0.2, -0.7);
/// assert!((res.re - -0.11504928).abs() < 1.0e-8 && (res.im - 0.09690468).abs() < 1.0e-8);
/// ```
pub fn spherical_harmonics(l: usize, m: i32, theta: f64, phi: f64) -> Complex {

    // We do the computation for the positive value
    let mp: i32 = m.abs();
    let cpx: Complex = Complex::from(0, mp as f64 * phi).exp();
    let poly = polynomial::Legendre::new(l, mp);

    // We follow QM norm
    let norm: f64 = (2 * l + 1) as f64 / (4.0 * PI);
    let top: f64 = basic::factorial(l - mp as usize) as f64;
    let bot: f64 = basic::factorial(l + mp as usize) as f64;

    // Computation with Legendre polynomial
    // (-1.0_f64).powi(m) term for the Condon-Shortley phase
    let res: Complex = (-1.0_f64).powi(m) * cpx * (norm * top / bot).sqrt() * poly.compute(theta.cos());

    // Modifying the result depending on the sign of m
    if m.is_negative() {
        res.conjugate() * (-1.0_f64).powi(mp)
    } else {
        res
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
