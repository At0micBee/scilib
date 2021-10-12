//!
//! # Spherical coordinates
//! 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

use std::fmt::{     // Formatter display
    Display,        // The display itself
    Result as DRes  // The associated result
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # Spherical coordinates
/// 
/// Defined for 3D space.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Spherical {
    /// Radial distance
    pub r: f64,
    /// Longitude angle
    pub theta: f64,
    /// Latitude angle
    pub phi: f64
}

/// # Display for Cartesian
/// 
/// Simply shows each value associated to an axis.
impl Display for Spherical {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> DRes {
        write!(f, "r={} :: theta={} :: phi={}", self.r, self.theta, self.phi)?;
        Ok(())
    }
}

impl Spherical {
    /// # Creates a new entity
    /// 
    /// Returns the same value as `Self::default()`, all elements are equal to zero.
    /// 
    /// ```
    /// # use scilib::coordinate::spherical::Spherical;
    /// let m = Spherical { r: 0.0, theta: 0.0, phi: 0.0 };
    /// let n = Spherical::new();
    /// let d = Spherical::default();
    /// 
    /// assert_eq!(m, n);
    /// assert_eq!(n, d);
    /// ```
    pub fn new() -> Self {
        Self::default()
    }

    /// # From the point
    /// 
    /// Creates a Cartesian struct from three given points in space.
    /// 
    /// ```
    /// # use scilib::coordinate::spherical::Spherical;
    /// let m = Spherical { r: 1.0, theta: 0.12, phi: 2.8 };
    /// let f = Spherical::from(1, 0.12, 2.8);
    /// 
    /// assert_eq!(m, f);
    /// ```
    pub fn from<T, U, V>(r: T, theta: U, phi: V) -> Self
    where T: Into<f64>, U: Into<f64>, V: Into<f64> {
        Self {
            r: r.into(),
            theta: theta.into(),
            phi: phi.into()
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
