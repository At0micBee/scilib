//!
//! # Methods for quantum mechanics
//! 

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

/// # Computing legal angular momentum numbers l based on n
/// 
/// This quantum number exists in the range [0, n-1].
/// 
/* /// ```
/// # use scilib::quantum::get_l;
/// let l = get_l(3);
/// assert_eq!(l, vec![0, 1, 2]);
/// ``` */
pub fn get_l<T>(n: T) -> Vec<usize>
where T: Into<usize> {
    (0..n.into()).collect()
}

/// # Computing legal magnetic quantum numbers m based on l
/// 
/// This quantum number exists in the range [-l, l].
/// 
/* /// ```
/// # use scilib::quantum::get_m;
/// let m = get_m(2);
/// assert_eq!(m, vec![-2, -1, 0, 1, 2]);
/// ``` */
pub fn get_m<T>(l: T) -> Vec<isize>
where T: Into<isize> {
    let li: isize = l.into();
    (-li..=li).collect()
}

/// # Radial wavefunctions
/// 
/// They are quite costly to derive properly, due to the normalization process. Luckily, we can 
/// compute them manually and implement them this way.
/// 
/// This function is a simple sub-function selector that return the result of the appropriate
/// radial function. The `match` allows great flexibility but incurs runtime penalty.
/// 
/// WARNING: Not all wavefunction have been implemented; max is n=2, l=1.
pub fn radial_wavefunction(n: usize, l: usize, r: f64) -> f64 {

    // Pre-computing
    let factor: f64 = r / (n as f64 * cst::A_0);

    // Computing the norm of the function
    let mut norm: f64 = (2.0 / (n as f64 * cst::A_0)).powi(3);
    norm *= basic::factorial(n - l - 1) as f64 / (2 * n * basic::factorial(n + l).pow(3)) as f64;

    // Computing the term associated to the Laguerre polynomial
    let poly = polynomial::Laguerre::new(n + 1, 2 * l as i32 + 1).compute(2.0 * factor);

    -(norm.sqrt()).powi(l as i32) * poly * (-factor.exp())
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # Spherical harmonics
/// 
/// Provides the solution for the Ylm spherical harmonics in quantum mechanics. Note the the Condon-Shortley
/// phase is computed withing the Legendre polynomial, this solution is also therefore viable for any acoustics.
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
    let res: Complex = cpx * (norm * top / bot).sqrt() * poly.compute(theta.cos());

    // Modifying the result depending on the sign of m
    if m.is_negative() {
        res.conjugate() * (-1.0_f64).powi(mp)
    } else {
        res
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
