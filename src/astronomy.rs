//!
//! # Methods for astronomy
//!

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

use std::f64::consts::{     // Using std lib constants
    PI                      // Pi
};

use std::fmt::{             // Formatter display
    Display,                // The display itself
    Result as DRes          // The associated result
};

use super::constant;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # Radec coordinate system
///
/// Right ascension and declination of the object in the sky. The values are stored as `f64` internally, and in radians.
/// The distance is kept in meters internally, for homogeneousness.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Radec {
    /// Right ascension
    pub ra: f64,
    /// Declination
    pub dec: f64,
    /// Distance from earth
    pub dist_earth: Option<f64>
}

/// # Display
///
/// Shows both ra and dec in degrees, which are more easily readable.
impl Display for Radec {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> DRes {
        write!(f, "ra={}° :: dec={}°", self.ra.to_degrees(), self.dec.to_degrees())?;

        // If there is a distance, we print it as well
        match self.dist_earth {
            Some(d) => write!(f, " :: distance={}ly", d / constant::LY)?,
            None => {}
        }
        Ok(())
    }
}

/// Implementing required methods
impl Radec {

    /// # Radec coordinates from angles in radians
    ///
    /// ```
    /// # use std::f64::consts::{ FRAC_PI_3, FRAC_PI_6 };
    /// # use scilib::astronomy::Radec;
    /// let coord1 = Radec::from_rad(FRAC_PI_6, FRAC_PI_3);
    /// let coord2 = Radec { ra: FRAC_PI_6, dec: FRAC_PI_3, dist_earth: None };
    /// assert!((coord1.ra - coord2.ra).abs() < 1.0e-15);
    /// assert!((coord1.dec - coord2.dec).abs() < 1.0e-15);
    /// ```
    pub fn from_rad<T, U>(ra: T, dec: U) -> Self
    where T: Into<f64>, U: Into<f64> {
        Self {
            ra: ra.into(),
            dec: dec.into(),
            dist_earth: None
        }
    }

    /// # Radec coordinates from angles in degrees
    ///
    /// ```
    /// # use std::f64::consts::{ FRAC_PI_3, FRAC_PI_6 };
    /// # use scilib::astronomy::Radec;
    /// let coord1 = Radec::from_degree(30.0, 60.0);
    /// let coord2 = Radec { ra: FRAC_PI_6, dec: FRAC_PI_3, dist_earth: None };
    /// assert!((coord1.ra - coord2.ra).abs() < 1.0e-15);
    /// assert!((coord1.dec - coord2.dec).abs() < 1.0e-15);
    /// ```
    pub fn from_degree<T, U>(ra: T, dec: U) -> Self
    where T: Into<f64>, U: Into<f64> {
        Self {
            ra: ra.into().to_radians(),
            dec: dec.into().to_radians(),
            dist_earth: None
        }
    }

    /// # Specifying known distance
    /// 
    /// Internal structure uses distance in meters.
    /// 
    /// ```
    /// # use scilib::constant;
    /// # use scilib::astronomy::Radec;
    /// let mut coord = Radec::from_degree(30.0, 12.3);
    /// coord.specify_distance(4.403 * constant::LY);
    /// assert_eq!(coord.ra, 30.0_f64.to_radians());
    /// assert_eq!(coord.dec, 12.3_f64.to_radians());
    /// assert_eq!(coord.dist_earth, Some(4.403 * constant::LY));
    /// ```
    pub fn specify_distance<T>(&mut self, distance: T)
    where T: Into<f64> {
        self.dist_earth = Some(distance.into());
    }

    /// # Angular separation computation
    ///
    /// We use the Vincenty formula to compute the angular separation, as it should remain stable no matter
    /// the separations and positions of the objects. The resulting separation will be in the range [0, pi];
    ///
    /// ```
    /// use scilib::astronomy::Radec;
    /// let c1 = Radec::from_rad(1.0, -1.2);
    /// let c2 = Radec::from_rad(-2.02, 0.13);
    ///
    /// assert!((c1.separation(c2) - 2.0685709648870154).abs() < 1.0e-15);
    /// ```
    pub fn separation(&self, other: Self) -> f64 {

        // Difference in right ascension
        let da: f64 = self.ra - other.ra;

        //Computing the two top terms
        let t1: f64 = (other.dec.cos() * da.sin()).powi(2);
        let t2: f64 = (self.dec.cos() * other.dec.sin() - self.dec.sin() * other.dec.cos() * da.cos()).powi(2);

        // Computing the bottom term
        let b: f64 = self.dec.sin() * other.dec.sin() + self.dec.cos() * other.dec.cos() * da.cos();

        let res: f64 = ((t1 + t2).sqrt() / b).atan();

        if res < 0.0 {
            res + PI
        } else {
            res
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # Apparent magnitude
/// Uses the standard reference luminosity to compute the apparent magnitude of an object.
/// $$
/// m_\mathrm{obj} = -2.5\log\left( \frac{b_\mathrm{obj}}{b_\mathrm{ref}} \right)
/// $$
/// Where $F_\mathrm{obj}$ is the flux of the object and $F_\mathrm{ref}$ is the reference flux (see constants).
/// 
/// ```
/// # use scilib::constant;
/// # use scilib::astronomy::apparent_mag;
/// let alpha_b = 0.4 * constant::SUN_L;        // Alpha centauri B, roughly
/// let distance = 4.403 * constant::LY;        // Seen from Earth
/// let res = apparent_mag(alpha_b, distance);
/// assert!((res - 1.36).abs() <= 0.03);
/// ```
pub fn apparent_mag(lum: f64, dist: f64) -> f64 {
    // We avoid a division by using the pre-computed value
    -2.5 * irradiance(lum, dist).log10() + constant::APP_MAG_SHIFT
}

/// # Absolute magnitude
/// Same approach that the apparent magnitude, but set at a distance of 10pc.
/// $$
/// M = -2.5\log\left( \frac{L}{L_0} \right)
/// $$
/// 
/// ```
/// # use scilib::constant;
/// # use scilib::astronomy::absolute_mag;
/// let abs = absolute_mag(constant::SUN_L);    // Sun absolute magnitude
/// assert!((abs - 4.74).abs() <= 0.02);        // As defined IAU
/// ```
pub fn absolute_mag(lum: f64) -> f64 {
    // We avoid a division by using the pre-computed value
    -2.5 * lum.log10() + constant::ABS_MAG_SHIFT
}

/// # Distance modulus
/// Computes the distance of an object based on its apparent and absolute magnitudes:
/// $$
/// d = 10^{1 + \frac{m-M}{5}}
/// $$
/// Where $m$ is the apparent magnitude and $M$ is the absolute magnitude.
/// 
/// ```
/// # use scilib::constant;
/// # use scilib::astronomy::{ apparent_mag, absolute_mag, distance_mod };
/// let app = apparent_mag(constant::SUN_L, 10.0 * constant::PC);
/// let abs = absolute_mag(constant::SUN_L);
/// let res = distance_mod(app, abs) / constant::PC;
/// assert!((res - 10.0).abs() <= 1e-5);
/// ```
pub fn distance_mod(m_app: f64, m_abs: f64) -> f64 {
    10.0_f64.powf(1.0 + (m_app - m_abs) / 5.0) * constant::PC
}

/// # Equilibrium temperature
///
/// Knowing the stellar effective temperature, its radius, the distance from the star
/// and the albedo of the considered planet, we can compute the equilibrium temperature
/// found at the planet.
pub fn t_eq(star_t: f64, star_rad: f64, dist: f64, albedo: f64) -> f64 {
    // Following the formula
    star_t * (star_rad / (2.0 * dist)).sqrt() * (1.0 - albedo).sqrt()
}

/// # Object irradiance
///
/// Computes the irradiance at a certain distance from an object, using its luminosity.
/// The `luminosity` is in Watts and the `distance` is in meters, result is in `W.m-2`.
///
/// ```rust
/// # use scilib::constant;
/// # use scilib::astronomy;
/// // Computing the solar irradiance at the Earth
/// let p_earth: f64 = astronomy::irradiance(constant::SUN_L, constant::AU);
/// assert!((p_earth - 1361.0).abs() < 1.0);
///
/// // Computing the irradiance from the sun at Alpha-Centauri
/// let p: f64 = astronomy::irradiance(constant::SUN_L, 4.3 * constant::LY);
/// assert!((p - 1.84e-8).abs() < 0.1e-8);
/// ```
pub fn irradiance(luminosity: f64, distance: f64) -> f64 {
    // Following the formula
    luminosity / (4.0 * PI * distance.powi(2))
}

/// # Planetary luminosity
///
/// Computes the luminosity of a planet, based on the received irradiance and the albedo.
pub fn planet_luminosity(albedo: f64, received: f64) -> f64 {
    received * albedo
}

/// # Energy received by an object
///
/// Computes the received energy by an object of a given surface, at a known distance.
/// Makes use of the `insolation` function to compute the surfacic power outputed by the object.
/// This function can be used both for primary sources of power, or reflective ones; the albedo
/// should be set to 1 for primary emitters.
pub fn received_energy(luminosity: f64, distance: f64, surface: f64) -> f64 {
    // Following the formula
    irradiance(luminosity, distance) * surface
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
