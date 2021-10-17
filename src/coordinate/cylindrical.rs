//!
//! # Cylindrical coordinates
//! 
//! Cylindrical coordinates use the distance to the origin of the point on a given plane, the angle to a reference on
//! that plane, and the altitude of the point.
//! - r: distance to origin, `[0, +∞`[
//! - theta: azimuth (longitude) of the point, `[0, 2π[`
//! - z: elevation (altitude) of the point, `]-∞, +∞`[
//! 
//! 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

use std::f64::consts::TAU;

use std::fmt::{     // Formatter display
    Display,        // The display itself
    Result as DRes  // The associated result
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # Cylindrical coordinates
/// 
/// Defined for 3D space
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Cylindrical {
    /// Radius from center of plane
    pub r: f64,
    /// Angle from reference on the plane
    pub theta: f64,
    /// Altitude
    pub z: f64
}


/// # Display for Cartesian
/// 
/// Simply shows each value associated to an axis.
impl Display for Cylindrical {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> DRes {
        write!(f, "r={} :: theta={}° :: z={}", self.r, self.theta.to_degrees(), self.z)?;
        Ok(())
    }
}

impl Cylindrical {
    /// # Creates a new entity
    /// 
    /// Returns the same value as `Self::default()`, all elements are equal to zero.
    /// 
    /// ```
    /// # use scilib::coordinate::cylindrical::Cylindrical;
    /// let m = Cylindrical { r: 0.0, theta: 0.0, z: 0.0 };
    /// let n = Cylindrical::new();
    /// let d = Cylindrical::default();
    /// 
    /// assert_eq!(m, n);
    /// assert_eq!(n, d);
    /// ```
    pub fn new() -> Self {
        Self::default()
    }

    /// # From the point
    /// 
    /// Creates a Cylindrical struct from three given points in space.
    /// 
    /// ```
    /// # use scilib::coordinate::cylindrical::Cylindrical;
    /// let m = Cylindrical { r: 1.0, theta: 0.12, z: 2.8 };
    /// let f = Cylindrical::from(1, 0.12, 2.8);
    /// 
    /// assert_eq!(m, f);
    /// ```
    pub fn from<T, U, V>(r: T, theta: U, z: V) -> Self
    where T: Into<f64>, U: Into<f64>, V: Into<f64> {
        Self {
            r: r.into(),
            theta: theta.into() % TAU,
            z: z.into()
        }
    }

    /// # From the point (degrees)
    /// 
    /// Creates a Cylindrical struct from three given points in space, with the angles in degrees.
    /// 
    /// ```
    /// # use scilib::coordinate::cylindrical::Cylindrical;
    /// let m = Cylindrical { r: 1.0, theta: 45.0_f64.to_radians(), z: -4.2 };
    /// let f = Cylindrical::from_degree(1, 45, -4.2);
    /// 
    /// assert_eq!(m, f);
    /// ```
    pub fn from_degree<T, U, V>(r: T, theta: U, z: V) -> Self
    where T: Into<f64>, U: Into<f64>, V: Into<f64> {
        let td: f64 = theta.into();
        Self {
            r: r.into(),
            theta: td.to_radians() % TAU,
            z: z.into()
        }
    }

    /// # Distance between two points
    /// 
    /// ```
    /// # use std::f64::consts::SQRT_2;
    /// # use scilib::coordinate::cylindrical::Cylindrical;
    /// let s1 = Cylindrical::from_degree(SQRT_2, 45, 0.0);
    /// let s2 = Cylindrical::from_degree(SQRT_2, -45, 0.0);
    /// 
    /// assert_eq!(s1.distance(s2), 2.0);
    /// ```
    pub fn distance(&self, other: Self) -> f64 {
        let t1: f64 = self.r.powi(2) + other.r.powi(2);
        let t2: f64 = (self.theta - other.theta).cos() * 2.0 * self.r * other.r;
        let t3: f64 = (self.z - other.z).powi(2);

        (t1 - t2 + t3).sqrt()
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
