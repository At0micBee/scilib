//!
//! # Methods for astronomy
//! 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # Radec coordinate system
/// 
/// Right ascension and declination of the object in the sky. The values are stored as `f64` internally, and in radians.
#[derive(Debug, Default, PartialEq)]
pub struct Radec {
    pub ra: f64,
    pub dec: f64
}

/// Implementing required methods
impl Radec {

    /// # Radec coordinates from angles in radians
    /// 
    /// ```
    /// # use std::f64::consts::{ FRAC_PI_3, FRAC_PI_6 };
    /// # use scilib::astronomy::Radec;
    /// let coord1 = Radec::from_rad(FRAC_PI_6, FRAC_PI_3);
    /// let coord2 = Radec { ra: FRAC_PI_6, dec: FRAC_PI_3 };
    /// assert!((coord1.ra - coord2.ra).abs() < 1.0e-15);
    /// assert!((coord1.dec - coord2.dec).abs() < 1.0e-15);
    /// ```
    pub fn from_rad(ra: f64, dec: f64) -> Self {
        Self {
            ra,
            dec
        }
    }

    /// # Radec coordinates from angles in degrees
    /// 
    /// ```
    /// # use std::f64::consts::{ FRAC_PI_3, FRAC_PI_6 };
    /// # use scilib::astronomy::Radec;
    /// let coord1 = Radec::from_degrees(30.0, 60.0);
    /// let coord2 = Radec { ra: FRAC_PI_6, dec: FRAC_PI_3 };
    /// assert!((coord1.ra - coord2.ra).abs() < 1.0e-15);
    /// assert!((coord1.dec - coord2.dec).abs() < 1.0e-15);
    /// ```
    pub fn from_degrees(ra: f64, dec: f64) -> Self {
        Self {
            ra: ra.to_radians(),
            dec: dec.to_radians()
        }
    }

    /// # Angular distance computation
    /// 
    /// We use the Vincenty formula to compute the angular distance, as it should remain stable no matter
    /// the distances and positions of the objects.
    pub fn distance(&self, other: Self) -> f64 {

        // Difference in right ascension
        let da: f64 = self.ra - other.ra;
        
        //Computing the two top terms
        let t1: f64 = (other.dec.cos() * da.sin()).powi(2);
        let t2: f64 = (self.dec.cos() * other.dec.sin() - self.dec.sin() * other.dec.cos() * da.cos()).powi(2);

        // Computing the bottom term
        let b: f64 = self.dec.sin() * other.dec.sin() + self.dec.cos() * other.dec.cos() * da.cos();
        
        ((t1 + t2).sqrt() / b).atan()
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
