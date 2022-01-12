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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # Radec coordinate system
///
/// Right ascension and declination of the object in the sky. The values are stored as `f64` internally, and in radians.
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
            Some(d) => write!(f, " :: distance={}ly", d)?,
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
