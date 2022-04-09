//!
//! # Planck laws & black body
//! 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

use super::constant as cst;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # Frequency Planck's law
/// ## Definition
/// This computes the following equation:
/// $$
/// B_\nu(\nu,T) = \frac{2h\nu^3}{c^2} \frac{1}{\exp\left(\frac{h\nu}{k_BT}\right) - 1}
/// $$
/// Where $h$ is the Planck constant, $c$ is the speed of light and $k_B$ is the Boltzmann constant.
/// 
/// ## Inputs
/// - `temperature`: temperature of the black body ($T$), in kelvin ($K$)
/// - `nu`: the frequency to compute ($\nu$), in meters ($m$)
/// 
/// Returns the corresponding irradiance.
/// 
/// ## Example
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
/// ## Definition
/// This computes the following equation:
/// $$
/// B_\nu(\nu,T) = \frac{2h\nu^3}{c^2} \frac{1}{\exp\left(\frac{h\nu}{k_BT}\right) - 1}
/// $$
/// Where $h$ is the Planck constant, $c$ is the speed of light and $k_B$ is the Boltzmann constant.
/// 
/// ## Inputs
/// - `temperature`: temperature of the black body ($T$), in kelvin ($K$)
/// - `nu`: the frequencies to compute ($\nu$), in meters ($m$)
/// 
/// Returns the corresponding irradiances.
/// 
/// ## Example
/// ```
/// # use scilib::range;
/// # use scilib::planck;
/// let freq = range::linear(0.0, 2000e12, 15);
/// let temp: f64 = 5700.0;
/// let res = planck::frequency_vec(temp, &freq);
/// assert!((res[8] - 1.458e-9).abs() < 1.0e-10);
/// ```
pub fn frequency_vec(temperature: f64, nu: &[f64]) -> Vec<f64> {

    let factor: f64 = 2.0 * cst::H / cst::C.powi(2);
    let pre: f64 = cst::H / (cst::K_B * temperature);

    nu.iter().map(|val| factor * val.powi(3) / ((pre * val).exp() - 1.0) ).collect()
}

/// # Wavelength Planck's law
/// ## Definition
/// This computes the following equation:
/// $$
/// B_\lambda(\lambda,T) = \frac{2hc^2}{\lambda^5} \frac{1}{\exp\left(\frac{hc}{\lambda k_BT}\right) - 1}
/// $$
/// Where $h$ is the Planck constant, $c$ is the speed of light and $k_B$ is the Boltzmann constant.
/// 
/// ## Inputs
/// - `temperature`: temperature of the black body ($T$), in kelvin ($K$)
/// - `lambda`: the wavelength to compute ($\lambda$), in meters ($m$)
/// 
/// Returns the corresponding irradiance.
/// 
/// ## Example
/// ```
/// # use scilib::planck;
/// let res = planck::wavelength(5700.0, 0.4e-6);
/// assert!((res - 2.118e13).abs() <= 1.0e10);
/// ```
pub fn wavelength(temperature: f64, lambda: f64) -> f64 {

    let factor: f64 = 2.0 * cst::H * cst::C.powi(2) / lambda.powi(5);

    factor / ((cst::H * cst::C / (lambda * cst::K_B * temperature)).exp() - 1.0)
}

/// # Wavelength Planck's law (spectrum)
/// ## Definition
/// This computes the following equation:
/// $$
/// B_\lambda(\lambda,T) = \frac{2hc^2}{\lambda^5} \frac{1}{\exp\left(\frac{hc}{\lambda k_BT}\right) - 1}
/// $$
/// Where $h$ is the Planck constant, $c$ is the speed of light and $k_B$ is the Boltzmann constant.
/// 
/// ## Inputs
/// - `temperature`: temperature of the black body ($T$), in kelvin ($K$)
/// - `lambda`: the wavelengths to compute ($\lambda$), in meters ($m$)
/// 
/// Returns the corresponding irradiances.
/// 
/// ```
/// # use scilib::range;
/// # use scilib::planck;
/// let wave = range::linear(0.0, 5e-6, 15);
/// let temp: f64 = 5700.0;
/// let res = planck::wavelength_vec(temp, &wave);
/// assert!((res[8] - 4.408e11).abs() < 1.0e8);
/// ```
pub fn wavelength_vec(temperature: f64, lambda: &[f64]) -> Vec<f64> {

    let factor: f64 = 2.0 * cst::H * cst::C.powi(2);
    let pre: f64 = cst::H * cst::C / (cst::K_B * temperature);

    lambda.iter().map(|val| factor / (val.powi(5) * ( (pre / val).exp() - 1.0 )) ).collect()
}

/// # Wavenumber Planck's law
/// ## Definition
/// This computes the following equation:
/// $$
/// B_{\tilde{n}}(\tilde{n},T) = 2hc^2\tilde{n}^3 \frac{1}{\exp\left(\frac{hc\tilde{n}}{k_BT}\right) - 1}
/// $$
/// Where $h$ is the Planck constant, $c$ is the speed of light and $k_B$ is the Boltzmann constant.
/// 
/// ## Inputs
/// - `temperature`: temperature of the black body ($T$), in kelvin ($K$)
/// - `number`: the wavenumber to compute ($\tilde{n}$), in 1/meters ($m^{-1}$)
/// 
/// Returns the corresponding irradiance.
/// 
/// ## Example
pub fn wavenumber(temperature: f64, number: f64) -> f64 {

    let factor: f64 = 2.0 * cst::H * cst::C.powi(2) * number.powi(3);

    factor / ((cst::H * cst::C * number / (cst::K_B * temperature)).exp() - 1.0)
}

/// # Wavenumber Planck's law (spectrum)
/// ## Definition
/// This computes the following equation:
/// $$
/// B_{\tilde{n}}(\tilde{n},T) = 2hc^2\tilde{n}^3 \frac{1}{\exp\left(\frac{hc\tilde{n}}{k_BT}\right) - 1}
/// $$
/// Where $h$ is the Planck constant, $c$ is the speed of light and $k_B$ is the Boltzmann constant.
/// 
/// ## Inputs
/// - `temperature`: temperature of the black body ($T$), in kelvin ($K$)
/// - `number`: the wavenumbers to compute ($\tilde{n}$), in 1/meters ($m^{-1}$)
/// 
/// Returns the corresponding irradiances.
/// 
/// ## Example
pub fn wavenumber_vec(temperature: f64, number: &[f64]) -> Vec<f64> {

    let factor: f64 = 2.0 * cst::H * cst::C.powi(2);
    let pre: f64 = cst::H * cst::C / (cst::K_B * temperature);

    number.iter().map(|val| factor * val.powi(3) / ( (pre * val).exp() - 1.0) ).collect()
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # Maximum emission of a black body (wavelength)
/// ## Definition
/// This is Wien's displacement law, applied to wavelength:
/// $$
/// \lambda_\mathrm{peak} = \frac{b}{T}
/// $$
/// Where $b$ is Wien's constant and $T$ is the temperature of the black body.
/// 
/// ## Inputs
/// - `temperature`: the temperature of the black body ($T$), in kelvin ($K$)
/// 
/// Returns the wavelength peak of the black body.
/// 
/// ## Example
/// ```
/// # use scilib::planck::peak_wave;
/// # use scilib::constant;
/// let sun_max = peak_wave(constant::SUN_TEFF);    // Max is greenish
/// assert!((sun_max - 5.0203949324324324e-7).abs() < 1.0e-16);
/// ```
pub fn peak_wave(temperature: f64) -> f64 {
    cst::WIEN_B / temperature
}

/// # Maximum emission of a black body (frequency)
/// ## Definition
/// This is Wien's displacement law, applied to frequency:
/// $$
/// \nu_\mathrm{peak} = b'T
/// $$
/// Where $b'$ is Wien's constant for frequency and $T$ is the temperature of the black body.
/// 
/// ## Inputs
/// - `temperature`: the temperature of the black body ($T$), in kelvin ($K$)
/// 
/// Returns the frequency peak of the black body.
/// 
/// ## Example
/// ```
/// # use scilib::planck::peak_freq;
/// # use scilib::constant;
/// let sun_max = peak_freq(constant::SUN_TEFF);    // Max is greenish
/// assert!((sun_max - 339331594694040.0).abs() < 1.0e-16);
/// ```
pub fn peak_freq(temperature: f64) -> f64 {
    cst::WIEN_B_FREQ * temperature
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
