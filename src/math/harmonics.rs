//!
//! Implementation of the Spherical harmonics function
//! 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

use std::f64::consts::PI;           // Value of pi

use super::{                        // Calling other modules
    basic,
    complex::Complex,
    polynomial
};   

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # Spherical harmonics
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
