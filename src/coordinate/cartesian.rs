//!
//! # Cartesian coordinates
//! 
//! Cartesian coordinates store the distance of the point compared to the origin for each axis.
//! 
//! Support conversion to and from Spherical and Cylindrical coordinates.

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

use std::f64::consts::{         // Using std lib constants
    PI                          // Pi
};

use std::ops::{                 // Implementing basic operations
    Add,                        // Addition
    AddAssign,                  // Assigning addition
    Sub,                        // Subtraction
    SubAssign,                  // Assigning addition
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
    cylindrical::Cylindrical,   // Cylindrical coordinates
    spherical::Spherical        // Spherical coordinates
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # Cartesian coordinates
/// 
/// Defined for 3D space.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Cartesian {
    /// x axis
    pub x: f64,
    /// y axis
    pub y: f64,
    /// z axis
    pub z: f64
}

/// # Display for Cartesian
/// 
/// Simply shows each value associated to an axis.
impl Display for Cartesian {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> DRes {
        write!(f, "x={} :: y={} :: z={}", self.x, self.y, self.z)?;
        Ok(())
    }
}

/// Implementing required methods
impl Cartesian {

    /// # Creates a new entity
    /// 
    /// Returns the same value as `Self::default()`, all elements are equal to zero.
    /// 
    /// ```
    /// # use scilib::coordinate::cartesian::Cartesian;
    /// let m = Cartesian { x: 0.0, y: 0.0, z: 0.0 };
    /// let n = Cartesian::origin();
    /// let d = Cartesian::default();
    /// 
    /// assert_eq!(m, n);
    /// assert_eq!(n, d);
    /// ```
    pub const fn origin() -> Self {
        Self {
            x: 0.0,
            y: 0.0,
            z: 0.0
        }
    }

    /// # From each point
    /// 
    /// Creates a Cartesian struct from three given points in space.
    /// 
    /// ```
    /// # use scilib::coordinate::cartesian::Cartesian;
    /// let m = Cartesian { x: 3.2, y: -2.0, z: 0.5e-3 };
    /// let f = Cartesian::from(3.2, -2, 0.5e-3);
    /// 
    /// assert_eq!(m, f);
    /// ```
    pub fn from<T, U, V>(x: T, y: U, z: V) -> Self
    where T: Into<f64>, U: Into<f64>, V: Into<f64> {
        Self {
            x: x.into(),
            y: y.into(),
            z: z.into()
        }
    }

    /// # From another coordinate system
    /// 
    /// Creates an Cartesian struct from another coordinate system. Calls the `Into<Cartesian>` method,
    /// which is verified in its implementation.
    /// 
    /// ```
    /// # use scilib::coordinate::cartesian::Cartesian;
    /// # use scilib::coordinate::spherical::Spherical;
    /// # use scilib::coordinate::cylindrical::Cylindrical;
    /// let s: Spherical = Spherical::from(2.4, 65, 15);
    /// let c: Cylindrical = Cylindrical::from_degree(1.2, 32, -3);
    /// let res1: Cartesian = Cartesian::from_coord(s);
    /// let res2: Cartesian = Cartesian::from_coord(c);
    /// 
    /// assert_eq!(res1, s.into());
    /// assert_eq!(res2, c.into());
    /// ```
    pub fn from_coord<T>(c: T) -> Self
    where T: Into<Self> {
        c.into()
    }

    /// # Computes the vector norm
    /// 
    /// We follow the convention of the l2 norm in this implementation.
    /// 
    /// ```
    /// # use scilib::coordinate::cartesian::Cartesian;
    /// let v = Cartesian::from(1, -2.5, 5.2);
    /// let n = v.norm();
    /// 
    /// assert!((n - 5.85576638878294).abs() < 1.0e-14);
    /// ```
    pub fn norm(&self) -> f64 {
        (self.x.powi(2) + self.y.powi(2) + self.z.powi(2)).sqrt()
    }

    /// # Distance squared between two points
    /// 
    /// Computes the square of the distance between two points. As the distance is the L2 norm,
    /// it requires taking the square root of the sum of the squared difference of each term.
    /// Computing the distance squared is therefore more efficient, as it simply skips computing
    /// the square root followed by the the squaring operation.
    /// 
    /// ```
    /// # use scilib::coordinate::cartesian::Cartesian;
    /// let f = Cartesian::from(1.0, 2, 0.5);
    /// let p = Cartesian::from(2, -1, 3.5);
    /// assert!((f.distance_square(p) - 19.0).abs() < 1.0e-15);
    /// ```
    pub fn distance_square(self, rhs: Self) -> f64 {
        (rhs.x - self.x).powi(2) + (rhs.y - self.y).powi(2) + (rhs.z - self.z).powi(2)
    }

    /// # Distance between two points
    /// 
    /// Computes the distance between two points in space.
    /// 
    /// ```
    /// # use scilib::coordinate::cartesian::Cartesian;
    /// let f = Cartesian::from(1.0, 2, 0.5);
    /// let p = Cartesian::from(2, -1, 3.5);
    /// assert!((f.distance(p) - 4.358898943540674).abs() < 1.0e-15);
    /// ```
    pub fn distance(self, rhs: Self) -> f64 {
        self.distance_square(rhs).sqrt()
    }

    /// # Coordinate rotation
    /// 
    /// Computes the resulting coordinates after an arbitrary rotation in 3D. The rotation
    /// arguments are the yaw, pitch and roll in radians.
    /// 
    /// - Yaw is the rotation around the `z` axis
    /// - Pitch is the rotation around the `y` axis
    /// - Roll is the rotation around the `x` axis
    /// 
    /// ```
    /// # use std::f64::consts::FRAC_PI_2;
    /// # use scilib::coordinate::cartesian::Cartesian;
    /// 
    /// // Setting a point in y=0, and rotating it 90° trigonometry-wise
    /// let f = Cartesian::from(0, 1, 0);
    /// let res = f.rotate(FRAC_PI_2, 0.0, 0.0);
    /// 
    /// assert_eq!(res.x, -1.0);
    /// assert!((res.y - 0.0).abs() < f64::EPSILON);
    /// assert_eq!(res.z, 0.0);
    /// ```
    pub fn rotate(&self, yaw: f64, pitch: f64, roll: f64) -> Self {

        let (y_sin, y_cos): (f64, f64) = yaw.sin_cos();
        let (p_sin, p_cos): (f64, f64) = pitch.sin_cos();
        let (r_sin, r_cos): (f64, f64) = roll.sin_cos();

        // First row of the matrix
        let a11: f64 = y_cos * p_cos;
        let a12: f64 = y_cos * p_sin * r_sin - y_sin * r_cos;
        let a13: f64 = y_cos * p_sin * r_cos + y_sin * r_sin;

        // Second row
        let a21: f64 = y_sin * p_cos;
        let a22: f64 = y_sin * p_sin * r_sin + y_cos * r_cos;
        let a23: f64 = y_sin * p_sin * r_cos - y_cos * r_sin;

        // Third row
        let a31: f64 = -p_sin;
        let a32: f64 = p_cos * r_sin;
        let a33: f64 = p_cos * r_cos;

        // Following matrix multiplication
        Self {
            x: a11 * self.x + a12 * self.y + a13 * self.z,
            y: a21 * self.x + a22 * self.y + a23 * self.z,
            z: a31 * self.x + a32 * self.y + a33 * self.z
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # Conversion to spherical coordinates
/// 
/// ```
/// # use std::f64::consts::SQRT_2;
/// # use scilib::coordinate::spherical::Spherical;
/// # use scilib::coordinate::cartesian::Cartesian;
/// 
/// let s = Cartesian::from(1, 1, 0);
/// let conv: Spherical = s.into();
/// let expected = Spherical::from_degree(SQRT_2, 45, 90);
/// 
/// assert_eq!(conv, expected);
/// 
/// let s = Cartesian::from(1, 0, 1);
/// let conv: Spherical = s.into();
/// let expected = Spherical::from_degree(SQRT_2, 0, 45);
/// 
/// assert_eq!(conv.r, expected.r);
/// assert_eq!(conv.theta, expected.theta);
/// assert!((conv.phi - expected.phi).abs() < 1.0e-15);
/// ```
impl Into<Spherical> for Cartesian {
    fn into(self) -> Spherical {
        let rho: f64 = self.norm();
        let mut nt: f64 = (self.y / self.x).atan();

        // If we were in the wrong quadrants, the atan range doesn't work
        if self.x.is_sign_negative() {
            nt += PI;
        }

        Spherical {
            r: rho,
            theta: nt,
            phi: (self.z / rho).acos()
        }
    }
}

/// # Conversion to cylindrical coordinates
/// 
/// ```
/// # use std::f64::consts::SQRT_2;
/// # use scilib::coordinate::cylindrical::Cylindrical;
/// # use scilib::coordinate::cartesian::Cartesian;
/// 
/// let s = Cartesian::from(1, 1, 1.2);
/// let conv: Cylindrical = s.into();
/// let expected = Cylindrical::from_degree(SQRT_2, 45, 1.2);
/// 
/// assert_eq!(conv, expected);
/// 
/// let s = Cartesian::from(1, 0, 1);
/// let conv: Cylindrical = s.into();
/// let expected = Cylindrical::from_degree(1.0, 0, 1.0);
/// 
/// assert_eq!(conv, expected);
/// ```
impl Into<Cylindrical> for Cartesian {
    fn into(self) -> Cylindrical {
        let rho: f64 = (self.x.powi(2) + self.y.powi(2)).sqrt();
        let mut nt: f64 = (self.y / self.x).atan();

        // If we were in the wrong quadrants, the atan range doesn't work
        if self.x.is_sign_negative() {
            nt += PI;
        }

        Cylindrical {
            r: rho,
            theta: nt,
            z: self.z
        }
    }
}

/// # Addition
/// 
/// ```
/// # use scilib::coordinate::cartesian::Cartesian;
/// let c1 = Cartesian::from(1.0, 2.5, 3.0);
/// let c2 = Cartesian::from(-1.0, 1, 0.0);
/// let res = c1 + c2;
/// let expected = Cartesian::from(0.0, 3.5, 3.0);
/// 
/// assert_eq!(res, expected);
/// ```
impl<T: Into<Self>> Add<T> for Cartesian {
    type Output = Self;
    fn add(self, rhs: T) -> Self::Output {
        let rhs: Self = rhs.into();
        Self {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z
        }
    }
}

/// # Addition with reference
/// 
/// ```
/// # use scilib::coordinate::cartesian::Cartesian;
/// let c1 = Cartesian::from(1.0, 2.5, 3.0);
/// let c2 = Cartesian::from(-1.0, 1, 0.0);
/// let res = c1 + &c2;
/// let expected = Cartesian::from(0.0, 3.5, 3.0);
/// 
/// assert_eq!(res, expected);
/// ```
impl Add<&Self> for Cartesian {
    type Output = Self;
    fn add(self, rhs: &Self) -> Self::Output {
        Self {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z
        }
    }
}

/// # Assigning addition
/// 
/// ```
/// # use scilib::coordinate::cartesian::Cartesian;
/// let mut c1 = Cartesian::from(1.0, 2.5, 3.0);
/// let c2 = Cartesian::from(-1.0, 1, 0.0);
/// c1 += c2;
/// let expected = Cartesian::from(0.0, 3.5, 3.0);
/// 
/// assert_eq!(c1, expected);
/// ```
impl<T: Into<Self>> AddAssign<T> for Cartesian {
    fn add_assign(&mut self, rhs: T) {
        let rhs: Self = rhs.into();
        self.x += rhs.x;
        self.y += rhs.y;
        self.z += rhs.z;
    }
}

/// # Assigning addition with reference
/// 
/// ```
/// # use scilib::coordinate::cartesian::Cartesian;
/// let mut c1 = Cartesian::from(1.0, 2.5, 3.0);
/// let c2 = Cartesian::from(-1.0, 1, 0.0);
/// c1 += &c2;
/// let expected = Cartesian::from(0.0, 3.5, 3.0);
/// 
/// assert_eq!(c1, expected);
/// ```
impl AddAssign<&Self> for Cartesian {
    fn add_assign(&mut self, rhs: &Self) {
        self.x += rhs.x;
        self.y += rhs.y;
        self.z += rhs.z;
    }
}

/// # Subtraction
/// 
/// ```
/// # use scilib::coordinate::cartesian::Cartesian;
/// let c1 = Cartesian::from(1.0, 2.5, 3.0);
/// let c2 = Cartesian::from(-1.0, 1, 0.0);
/// let res = c1 - c2;
/// let expected = Cartesian::from(2, 1.5, 3.0);
/// 
/// assert_eq!(res, expected);
/// ```
impl<T: Into<Self>> Sub<T> for Cartesian {
    type Output = Self;
    fn sub(self, rhs: T) -> Self::Output {
        let rhs: Self = rhs.into();
        Self {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z
        }
    }
}

/// # Subtraction with reference
/// 
/// ```
/// # use scilib::coordinate::cartesian::Cartesian;
/// let c1 = Cartesian::from(1.0, 2.5, 3.0);
/// let c2 = Cartesian::from(-1.0, 1, 0.0);
/// let res = c1 - &c2;
/// let expected = Cartesian::from(2, 1.5, 3.0);
/// 
/// assert_eq!(res, expected);
/// ```
impl Sub<&Self> for Cartesian {
    type Output = Self;
    fn sub(self, rhs: &Self) -> Self::Output {
        Self {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z
        }
    }
}

/// # Assigning subtraction
/// 
/// ```
/// # use scilib::coordinate::cartesian::Cartesian;
/// let mut c1 = Cartesian::from(1.0, 2.5, 3.0);
/// let c2 = Cartesian::from(-1.0, 1, 0.0);
/// c1 -= c2;
/// let expected = Cartesian::from(2, 1.5, 3.0);
/// 
/// assert_eq!(c1, expected);
impl<T: Into<Self>> SubAssign<T> for Cartesian {
    fn sub_assign(&mut self, rhs: T) {
        let rhs: Self = rhs.into();
        self.x -= rhs.x;
        self.y -= rhs.y;
        self.z -= rhs.z;
    }
}

/// # Assigning subtraction with reference
/// 
/// ```
/// # use scilib::coordinate::cartesian::Cartesian;
/// let mut c1 = Cartesian::from(1.0, 2.5, 3.0);
/// let c2 = Cartesian::from(-1.0, 1, 0.0);
/// c1 -= &c2;
/// let expected = Cartesian::from(2, 1.5, 3.0);
/// 
/// assert_eq!(c1, expected);
impl SubAssign<&Self> for Cartesian {
    fn sub_assign(&mut self, rhs: &Self) {
        self.x -= rhs.x;
        self.y -= rhs.y;
        self.z -= rhs.z;
    }
}

/// # Scalar multiplication
/// 
/// Multiplies each field by the scalar.
/// 
/// ```
/// # use scilib::coordinate::cartesian::Cartesian;
/// let c1 = Cartesian::from(-1.0, 2.5, 3.0);
/// let res = c1 * 2;
/// let expected = Cartesian::from(-2, 5, 6.0);
/// 
/// assert_eq!(res, expected);
/// ```
impl<T: Into<f64>> Mul<T> for Cartesian {
    type Output = Self;
    fn mul(self, rhs: T) -> Self::Output {
        let f: f64 = rhs.into();
        Self {
            x: self.x * f,
            y: self.y * f,
            z: self.z * f
        }
    }
}

/// # Assigning scalar multiplication
/// 
/// Multiplies each field by the scalar in place.
/// 
/// ```
/// # use scilib::coordinate::cartesian::Cartesian;
/// let mut c1 = Cartesian::from(-1.0, 2.5, 3.0);
/// c1 *= 2;
/// let expected = Cartesian::from(-2, 5, 6.0);
/// 
/// assert_eq!(c1, expected);
/// ```
impl<T: Into<f64>> MulAssign<T> for Cartesian {
    fn mul_assign(&mut self, rhs: T) {
        let f: f64 = rhs.into();
        self.x *= f;
        self.y *= f;
        self.z *= f;
    }
}

/// # Scalar division
/// 
/// Divises each field by the scalar.
/// 
/// ```
/// # use scilib::coordinate::cartesian::Cartesian;
/// let c1 = Cartesian::from(-2, 5, 6.0);
/// let res = c1 / 2;
/// let expected = Cartesian::from(-1.0, 2.5, 3.0);
/// 
/// assert_eq!(res, expected);
/// ```
impl<T: Into<f64>> Div<T> for Cartesian {
    type Output = Self;
    fn div(self, rhs: T) -> Self::Output {
        let f: f64 = 1.0 / rhs.into();
        Self {
            x: self.x * f,
            y: self.y * f,
            z: self.z * f,
        }
    }
}

/// # Assigning scalar division
/// 
/// Divises each field by the scalar in place.
/// 
/// ```
/// # use scilib::coordinate::cartesian::Cartesian;
/// let mut c1 = Cartesian::from(-2, 5, 6.0);
/// c1 /= 2;
/// let expected = Cartesian::from(-1.0, 2.5, 3.0);
/// 
/// assert_eq!(c1, expected);
/// ```
impl<T: Into<f64>> DivAssign<T> for Cartesian {
    fn div_assign(&mut self, rhs: T) {
        let f: f64 = 1.0 / rhs.into();
        self.x *= f;
        self.y *= f;
        self.z *= f;
    }
}

/// # Negation
/// 
/// Changes the values to their opposites.
/// 
/// ```
/// # use scilib::coordinate::cartesian::Cartesian;
/// let c1 = Cartesian::from(-2, 5, 6.0);
/// let c2 = -c1;
/// let expected = Cartesian::from(2, -5.0, -6);
/// 
/// assert_eq!(c2, expected);
/// ```
impl Neg for Cartesian {
    type Output = Self;
    fn neg(self) -> Self::Output {
        self * -1
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
