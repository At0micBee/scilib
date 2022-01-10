//!
//! # Methods for astronomy
//!

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

use std::f64::consts::PI;   // Calling constants

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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
