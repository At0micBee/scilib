//!
//! # Methods for quantum mechanics
//! 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

use std::f64::consts::PI;           // Value of pi

use crate::{                        // Calling other modules
    math::{                  
        basic,
        complex::Complex,
        polynomial
    },
    constant as cst
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # Computing legal angular momentum numbers l based on n
/// 
/// This quantum number exists in the range [0, n-1].
/// 
/// ```
/// # use scilib::quantum::get_l;
/// let l = get_l(3);
/// assert_eq!(l, vec![0, 1, 2]);
/// ```
#[allow(unused)]
fn get_l(n: usize) -> Vec<usize> {
    (0..n).collect()
}

/// # Computing legal magnetic quantum numbers m based on l
/// 
/// This quantum number exists in the range [-l, l].
/// 
/// ```
/// # use scilib::quantum::get_m;
/// let m = get_m(2);
/// assert_eq!(m, vec![-2, -1, 0, 1, 2]);
/// ```
#[allow(unused)]
fn get_m(l: usize) -> Vec<isize> {
    let li: isize = l as isize;
    (-li..=li).collect()
}

/// # Radial wavefunctions
/// 
/// They are quite costly to derive properly, due to the normalization process. Luckily, we can 
/// compute them manually and implement them this way.
#[allow(unused)]
fn radial_wavefunction(n: usize, l: usize, r: f64) -> f64 {

    match (n, l) {
        // n=1, l=0
        (1, 0) => {
            2.0 * (1.0 / cst::A_0).powf(3.0 / 2.0) * (-r / cst::A_0).exp()
        },
        // n=2, l=0
        (2, 0) => {
            2.0 * (0.5 / cst::A_0).powf(3.0 / 2.0) * (1.0 - r / (2.0 * cst::A_0)) * (-r / (2.0 * cst::A_0)).exp()
        },
        // n=2, l=1
        (2, 1) => {
            (1.0 / 3.0_f64.sqrt()) * (0.5 / cst::A_0).powf(3.0 / 2.0) * (r / cst::A_0) * (-r / (2.0 * cst::A_0)).exp()
        }
        // Either impossible case, or not computed
        _ => 0.0
    }
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
