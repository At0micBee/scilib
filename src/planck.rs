//!
//! # Planck laws & black body
//! 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

use super::constant as cst;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # Frequency Planck's law
pub fn frequency(temperature: f64, nu: f64) -> f64 {

    // 2 h n³ / c²
    let factor: f64 = 2.0 * nu.powi(3) * cst::H / cst::C.powi(2);

    // 1 / (exp(h n / k T) - 1)
    factor / (cst::H * nu / (cst::K_B * temperature)).exp() - 1.0
}

/// # Wavelength Planck's law
pub fn wavelength(temperature: f64, lambda: f64) -> f64 {

    // 2 h c² / l⁵
    let factor: f64 = 2.0 * cst::H * cst::C.powi(2) / lambda.powi(5);

    // 1 / (exp(h c / l k T) - 1)
    factor / (cst::H * cst::C / (lambda * cst::K_B * temperature)).exp() - 1.0
}

/// # Wavenumber Planck's law
pub fn wavenumber(temperature: f64, number: f64) -> f64 {

    // 2 h c² n³
    let factor: f64 = 2.0 * cst::H * cst::C.powi(2) * number.powi(3);

    // 1 / (exp(h c n / k T) - 1)
    factor / (cst::H * cst::C * number/ (cst::K_B * temperature)).exp() - 1.0
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
