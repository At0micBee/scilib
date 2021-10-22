//!
//! # Cylindrical coordinates
//! 
//! Cylindrical coordinates use the distance to the origin of the point on a given plane, the angle to a reference on
//! that plane, and the altitude of the point.
//! - r: distance to origin, `[0, +∞`[
//! - theta: azimuth (longitude) of the point, `[0, 2π[`
//! - z: elevation (altitude) of the point, `]-∞, +∞`[
//! 
//! Support conversion to and from Cartesian and Spherical coordinates.

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

use std::f64::consts::{         // Using std lib constants
    PI,                         // Pi
    TAU                         // Tau
};

use std::ops::{                 // Implementing basic operations
    Add,                        // Addition
    Sub,                        // Subtraction
    Mul,                        // Multiplication
    MulAssign,                  // Assigning multiplication
    Div,                        // Division
    DivAssign,                  // Assigning division
    Neg                         // Negation
};

use std::fmt::{                 // Formatter display
    Display,                    // The display itself
    Result as DRes              // The associated result
};

use super::{                    // Using parts from the crate
    cartesian::Cartesian,       // Cartesian coordinates
    spherical::Spherical        // Spherical coordinates
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

    /// # From another coordinate system
    /// 
    /// Creates an Cylindrical struct from another coordinate system. Calls the `Into<Cylindrical>` method,
    /// which is verified in its implementation.
    /// 
    /// ```
    /// # use scilib::coordinate::cartesian::Cartesian;
    /// # use scilib::coordinate::spherical::Spherical;
    /// # use scilib::coordinate::cylindrical::Cylindrical;
    /// let c: Cartesian = Cartesian::from(0, 12, 3.2);
    /// let s: Spherical = Spherical::from_degree(1.2, 32, 60);
    /// let res1: Cylindrical = Cylindrical::from_coord(c);
    /// let res2: Cylindrical = Cylindrical::from_coord(s);
    /// 
    /// assert_eq!(res1, c.into());
    /// assert_eq!(res2, s.into());
    /// ```
    pub fn from_coord<T>(c: T) -> Self
    where T: Into<Self> {
        c.into()
    }

    /// # Distance between two points
    /// 
    /// ```
    /// # use std::f64::consts::SQRT_2;
    /// # use scilib::coordinate::cylindrical::Cylindrical;
    /// let s1 = Cylindrical::from_degree(SQRT_2, 45, 1.0);
    /// let s2 = Cylindrical::from_degree(SQRT_2, -45, -1.0);
    /// 
    /// assert_eq!(s1.distance(s2), 2.0 * SQRT_2);
    /// ```
    pub fn distance(&self, other: Self) -> f64 {
        let t1: f64 = self.r.powi(2) + other.r.powi(2);
        let t2: f64 = (self.theta - other.theta).cos() * 2.0 * self.r * other.r;
        let t3: f64 = (self.z - other.z).powi(2);

        (t1 - t2 + t3).sqrt()
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # Conversion to cartesian coordinates
/// 
/// ```
/// # use scilib::coordinate::cartesian::Cartesian;
/// # use scilib::coordinate::cylindrical::Cylindrical;
/// let s = Cylindrical::from_degree(3, 60, 4);
/// let conv: Cartesian = s.into();
/// let expected = Cartesian::from(1.5, 2.598076211, 4);
/// 
/// assert!((conv.x - expected.x).abs() < 1.0e-9);
/// assert!((conv.y - expected.y).abs() < 1.0e-9);
/// assert_eq!(conv.z, expected.z);
/// 
/// let s = Cylindrical::from_degree(5, 30, -2);
/// let conv: Cartesian = s.into();
/// let expected = Cartesian::from(4.330127019, 2.5, -2);
/// 
/// assert!((conv.x - expected.x).abs() < 1.0e-9);
/// assert!((conv.y - expected.y).abs() < 1.0e-9);
/// assert_eq!(conv.z, expected.z);
/// ```
impl Into<Cartesian> for Cylindrical {
    fn into(self) -> Cartesian {
        Cartesian {
            x: self.r * self.theta.cos(),
            y: self.r * self.theta.sin(),
            z: self.z
        }
    }
}

/// # Conversion to spherical coordinates
/// 
/// ```
/// # use scilib::coordinate::spherical::Spherical;
/// # use scilib::coordinate::cylindrical::Cylindrical;
/// 
/// let s = Cylindrical::from_degree(3, 60, 4);
/// let conv: Spherical = s.into();
/// let expected = Spherical::from_degree(5, 60, 36.86989765);
/// 
/// assert_eq!(conv.r, expected.r);
/// assert_eq!(conv.theta, expected.theta);
/// assert!((conv.phi - expected.phi).abs() < 1.0e-9);
/// 
/// let s = Cylindrical::from_degree(3, 60, -4);
/// let conv: Spherical = s.into();
/// let expected = Spherical::from_degree(5, 60, 143.1301024);
/// 
/// assert_eq!(conv.r, expected.r);
/// assert_eq!(conv.theta, expected.theta);
/// assert!((conv.phi - expected.phi).abs() < 1.0e-9);
/// ```
impl Into<Spherical> for Cylindrical {
    fn into(self) -> Spherical {
        let rho: f64 = (self.r.powi(2) + self.z.powi(2)).sqrt();
        let mut np: f64 = (self.r / self.z).atan();

        if self.z.is_sign_negative() {
            np += PI;
        }

        Spherical {
            r: rho,
            theta: self.theta,
            phi: np
        }
    }
}

/// # Addition
/// 
/// Converts the coordinate in cartesian for addition, then returns them as Cylindrical.
impl<T: Into<Cartesian>> Add<T> for Cylindrical {
    type Output = Self;
    fn add(self, rhs: T) -> Self::Output {
        let s: Cartesian = self.into();
        let r: Cartesian = rhs.into();
        (s + r).into()
    }
}

/// # Subtraction
/// 
/// Converts the coordinate in cartesian for subtraction, then returns them as Cylindrical.
impl<T: Into<Cartesian>> Sub<T> for Cylindrical {
    type Output = Self;
    fn sub(self, rhs: T) -> Self::Output {
        let s: Cartesian = self.into();
        let r: Cartesian = rhs.into();
        (s - r).into()
    }
}

/// # Scalar multiplication
/// 
/// Multiplies the radius by a scalar.
/// 
/// ```
/// # use scilib::coordinate::cylindrical::Cylindrical;
/// let s = Cylindrical::from_degree(1, 180, 2.1);
/// let res = s * -2;
/// let expected = Cylindrical::from(2, 0, -4.2);
/// 
/// assert_eq!(res, expected);
/// ```
impl<T: Into<f64>> Mul<T> for Cylindrical {
    type Output = Self;
    fn mul(self, rhs: T) -> Self::Output {
        let f: f64 = rhs.into();
        let res: Self = Self {
            r: self.r * f.abs(),
            theta: self.theta,
            z: self.z * f.abs()
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
/// # use scilib::coordinate::cylindrical::Cylindrical;
/// let mut s = Cylindrical::from_degree(1, 180, 2.1);
/// s *= -2;
/// let expected = Cylindrical::from_degree(2, 0, -4.2);
/// 
/// assert_eq!(s, expected);
/// ```
impl<T: Into<f64>> MulAssign<T> for Cylindrical {
    fn mul_assign(&mut self, rhs: T) {
        let f: f64 = rhs.into();
        self.r *= f.abs();
        self.z *= f;

        // If the sign is negative, we need to flip the angles
        if f.is_sign_negative() {
            self.theta = (self.theta + PI) % TAU;
        }
    }
}

/// # Scalar division
/// 
/// Divides the radius by a scalar.
/// 
/// ```
/// # use scilib::coordinate::cylindrical::Cylindrical;
/// let s = Cylindrical::from_degree(2, 60, 2.1);
/// let res = s / 2;
/// let expected = Cylindrical::from_degree(1, 60, 1.05);
/// 
/// assert_eq!(res, expected);
/// ```
impl<T: Into<f64>> Div<T> for Cylindrical {
    type Output = Self;
    fn div(self, rhs: T) -> Self::Output {
        let f: f64 = rhs.into();
        let res: Self = Self {
            r: self.r / f.abs(),
            theta: self.theta,
            z: self.z / f.abs()
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
/// # use scilib::coordinate::cylindrical::Cylindrical;
/// let mut s = Cylindrical::from_degree(2, 30, 2.1);
/// s /= 2;
/// let expected = Cylindrical::from_degree(1, 30, 1.05);
/// 
/// assert_eq!(s, expected);
/// ```
impl<T: Into<f64>> DivAssign<T> for Cylindrical {
    fn div_assign(&mut self, rhs: T) {
        let f: f64 = rhs.into();
        self.r /= f.abs();
        self.z /= f;

        // If the sign is negative, we need to flip the angles
        if f.is_sign_negative() {
            self.theta = (self.theta + PI) % TAU;
        }
    }
}

/// # Negation
/// 
/// Going to the opposite point.
/// 
/// ```
/// # use scilib::coordinate::cylindrical::Cylindrical;
/// let c1 = Cylindrical::from_degree(2, 35, 6.0);
/// let c2 = -c1;
/// let expected = Cylindrical::from_degree(2, 215, -6);
/// 
/// assert_eq!(c2, expected);
/// ```
impl Neg for Cylindrical {
    type Output = Self;
    fn neg(self) -> Self::Output {
        Self {
            r: self.r,
            theta: (self.theta + PI) % TAU,
            z: -self.z
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
