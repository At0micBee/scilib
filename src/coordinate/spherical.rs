//!
//! # Spherical coordinates
//! 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

use std::f64::consts::{
    PI,
    FRAC_PI_2
};

use std::ops::{     // Implementing basic operations
    Add,            // Addition
    AddAssign,      // Assigning addition
    Sub,            // Subtraction
    SubAssign,      // Assigning addition
    Mul,            // Multiplication
    MulAssign,      // Assigning multiplication
    Div,            // Division
    DivAssign,      // Assigning division
    Neg             // Negation
};

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

/// # Scalar multiplication
/// 
/// Multiplies the radius by a scalar.
/// 
/// ```
/// # use scilib::coordinate::spherical::Spherical;
/// let s = Spherical::from(1, 0.2, 2.1);
/// let res = s * 2;
/// let expected = Spherical::from(2, 0.2, 2.1);
/// 
/// assert_eq!(res, expected);
/// ```
impl<T: Into<f64>> Mul<T> for Spherical {
    type Output = Self;
    fn mul(self, rhs: T) -> Self::Output {
        let f: f64 = rhs.into();
        let res: Self = Self {
            r: self.r * f.abs(),
            theta: self.theta,
            phi: self.phi
        };

        // If the sign is negative, we need to flip the angles
        if f.is_sign_negative() {
            -res
        } else {
            res
        }
    }
}

/// # Assigning scalar multiplication
/// 
/// Multiplies the radius by a scalar in place.
/// 
/// ```
/// # use scilib::coordinate::spherical::Spherical;
/// let mut s = Spherical::from(1, 0.2, 2.1);
/// s *= 2;
/// let expected = Spherical::from(2, 0.2, 2.1);
/// 
/// assert_eq!(s, expected);
/// ```
impl<T: Into<f64>> MulAssign<T> for Spherical {
    fn mul_assign(&mut self, rhs: T) {
        let f: f64 = rhs.into();
        self.r *= f.abs();

        // If the sign is negative, we need to flip the angles
        if f.is_sign_negative() {
            self.theta += PI;
            self.phi += FRAC_PI_2;
        }
    }
}

/// # Scalar division
/// 
/// Divides the radius by a scalar.
/// 
/// ```
/// # use scilib::coordinate::spherical::Spherical;
/// let s = Spherical::from(2, 0.2, 2.1);
/// let res = s / 2;
/// let expected = Spherical::from(1, 0.2, 2.1);
/// 
/// assert_eq!(res, expected);
/// ```
impl<T: Into<f64>> Div<T> for Spherical {
    type Output = Self;
    fn div(self, rhs: T) -> Self::Output {
        let f: f64 = rhs.into();
        let res: Self = Self {
            r: self.r / f.abs(),
            theta: self.theta,
            phi: self.phi
        };

        // If the sign is negative, we need to flip the angles
        if f.is_sign_negative() {
            -res
        } else {
            res
        }
    }
}

/// # Assigning scalar division
/// 
/// Divides the radius by a scalar in place.
/// 
/// ```
/// # use scilib::coordinate::spherical::Spherical;
/// let mut s = Spherical::from(2, 0.2, 2.1);
/// s /= 2;
/// let expected = Spherical::from(1, 0.2, 2.1);
/// 
/// assert_eq!(s, expected);
/// ```
impl<T: Into<f64>> DivAssign<T> for Spherical {
    fn div_assign(&mut self, rhs: T) {
        let f: f64 = rhs.into();
        self.r /= f.abs();

        // If the sign is negative, we need to flip the angles
        if f.is_sign_negative() {
            self.theta += PI;
            self.phi += FRAC_PI_2;
        }
    }
}

/// # Negation
/// 
/// Going to the opposite point.
/// 
/// ```
/// # use scilib::coordinate::cartesian::Cartesian;
/// let c1 = Cartesian::from(-2, 5, 6.0);
/// let c2 = -c1;
/// let expected = Cartesian::from(2, -5.0, -6);
/// 
/// assert_eq!(c2, expected);
/// ```
impl Neg for Spherical {
    type Output = Self;
    fn neg(self) -> Self::Output {
        Self {
            r: self.r,
            theta: self.theta + PI,
            phi: self.phi + FRAC_PI_2
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
