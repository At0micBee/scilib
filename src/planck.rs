//!
//! # Planck laws & black body
//! 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

use super::constant as cst;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # Frequency Planck's law
/// This computes the following equation:
/// $$
/// B_\nu(\nu,T) = \frac{2h\nu^3}{c^2} \frac{1}{\exp\left(\frac{h\nu}{k_BT}\right) - 1}
/// $$
/// Where $h$ is the Planck constant, $c$ is the speed of light and $k_B$ is the Boltzmann constant.
/// 
/// The function uses `temperature` ($T$) the temperature of the body, `nu` ($\nu$) the frequency to evaluate,
/// and returns the corresponding irradiance.
/// 
/// ```
/// # use scilib::planck;
/// let res = planck::frequency(5700.0, 650e12);
/// assert!((res - 1.708e-8).abs() < 1.0e-10);
/// ```
pub fn frequency(temperature: f64, nu: f64) -> f64 {

    let factor: f64 = 2.0 * nu.powi(3) * cst::H / cst::C.powi(2);

    factor / ((cst::H * nu / (cst::K_B * temperature)).exp() - 1.0)
}

/// # Frequency Planck's law (spectrum)
/// This computes the following equation:
/// $$
/// B_\nu(\nu,T) = \frac{2h\nu^3}{c^2} \frac{1}{\exp\left(\frac{h\nu}{k_BT}\right) - 1}
/// $$
/// Where $h$ is the Planck constant, $c$ is the speed of light and $k_B$ is the Boltzmann constant.
/// 
/// The function uses `temperature` ($T$) the temperature of the body, `nu` ($\nu$) the frequencies to evaluate,
/// and returns the corresponding irradiances.
/// 
/// ```
/// # use scilib::range;
/// # use scilib::planck;
/// let freq = range::linear(0.0, 2000e12, 15);
/// let temp: f64 = 5700.0;
/// let res = planck::frequency_vec(temp, &freq);
/// assert!((res[8] - 1.458e-9).abs() < 1.0e-10);
/// ```
pub fn frequency_vec(temperature: f64, nu: &Vec<f64>) -> Vec<f64> {

    let factor: f64 = 2.0 * cst::H / cst::C.powi(2);
    let pre: f64 = cst::H / (cst::K_B * temperature);

    nu.iter().map(|val| factor * val.powi(3) / ((pre * val).exp() - 1.0) ).collect()
}

/// # Wavelength Planck's law
/// This computes the following equation:
/// $$
/// B_\lambda(\lambda,T) = \frac{2hc^2}{\lambda^5} \frac{1}{\exp\left(\frac{hc}{\lambda k_BT}\right) - 1}
/// $$
/// Where $h$ is the Planck constant, $c$ is the speed of light and $k_B$ is the Boltzmann constant.
/// 
/// The function uses `temperature` ($T$) the temperature of the body, `lambda` ($\lambda$) the wavelength to evaluate,
/// and returns the corresponding irradiance.
pub fn wavelength(temperature: f64, lambda: f64) -> f64 {

    let factor: f64 = 2.0 * cst::H * cst::C.powi(2) / lambda.powi(5);

    factor / ((cst::H * cst::C / (lambda * cst::K_B * temperature)).exp() - 1.0)
}

/// # Wavelength Planck's law (spectrum)
/// This computes the following equation:
/// $$
/// B_\lambda(\lambda,T) = \frac{2hc^2}{\lambda^5} \frac{1}{\exp\left(\frac{hc}{\lambda k_BT}\right) - 1}
/// $$
/// Where $h$ is the Planck constant, $c$ is the speed of light and $k_B$ is the Boltzmann constant.
/// 
/// The function uses `temperature` ($T$) the temperature of the body, `lambda` ($\lambda$) the wavelengths to evaluate,
/// and returns the corresponding irradiances.
pub fn wavelength_vec(temperature: f64, lambda: &Vec<f64>) -> Vec<f64> {

    let factor: f64 = 2.0 * cst::H * cst::C.powi(2);
    println!("{}", factor);
    let pre: f64 = cst::H * cst::C / (cst::K_B * temperature);
    println!("{}", pre);

    lambda.iter().map(|val| factor / (val.powi(5) * ( (pre / val).exp() - 1.0 )) ).collect()
}

/// # Wavenumber Planck's law
/// This computes the following equation:
/// $$
/// B_{\tilde{n}}(\tilde{n},T) = 2hc^2\tilde{n}^3 \frac{1}{\exp\left(\frac{hc\tilde{n}}{k_BT}\right) - 1}
/// $$
/// Where $h$ is the Planck constant, $c$ is the speed of light and $k_B$ is the Boltzmann constant.
/// 
/// The function uses `temperature` ($T$) the temperature of the body, `number` ($\tilde{n}$) the wavenumber to evaluate,
/// and returns the corresponding irradiance.
pub fn wavenumber(temperature: f64, number: f64) -> f64 {

    let factor: f64 = 2.0 * cst::H * cst::C.powi(2) * number.powi(3);

    factor / ((cst::H * cst::C * number / (cst::K_B * temperature)).exp() - 1.0)
}

/// # Wavenumber Planck's law (spectrum)
/// This computes the following equation:
/// $$
/// B_{\tilde{n}}(\tilde{n},T) = 2hc^2\tilde{n}^3 \frac{1}{\exp\left(\frac{hc\tilde{n}}{k_BT}\right) - 1}
/// $$
/// Where $h$ is the Planck constant, $c$ is the speed of light and $k_B$ is the Boltzmann constant.
/// 
/// The function uses `temperature` ($T$) the temperature of the body, `number` ($\tilde{n}$) the wavenumbers to evaluate,
/// and returns the corresponding irradiances.
pub fn wavenumber_vec(temperature: f64, number: &Vec<f64>) -> Vec<f64> {

    let factor: f64 = 2.0 * cst::H * cst::C.powi(2);
    let pre: f64 = cst::H * cst::C / (cst::K_B * temperature);

    number.iter().map(|val| factor * val.powi(3) / ( (pre * val).exp() - 1.0) ).collect()
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
