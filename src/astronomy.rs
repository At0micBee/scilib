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
    pub ra: f64,
    pub dec: f64
}

/// # Display
/// 
/// Shows both ra and dec in degrees, which are more easily readable.
impl Display for Radec {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> DRes {
        write!(f, "ra={}° :: dec={}°", self.ra.to_degrees(), self.dec.to_degrees())?;
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
    /// the distances and positions of the objects. The resulting distance will be in the range [0, pi];
    /// 
    /// ```
    /// use scilib::astronomy::Radec;
    /// let c1 = Radec::from_rad(1.0, -1.2);
    /// let c2 = Radec::from_rad(-2.02, 0.13);
    /// 
    /// assert!((c1.distance(c2) - 2.0685709648870154).abs() < 1.0e-15);
    /// ```
    pub fn distance(&self, other: Self) -> f64 {

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
