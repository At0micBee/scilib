//!
//! # Coordinates systems
//! 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # Cartesian coordinates
/// 
/// Defined for 3D space.
#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub struct Cartesian {
    pub x: f64,
    pub y: f64,
    pub z: f64
}

impl Cartesian {

    /// # Creates a new entity
    /// 
    /// Returns the same value as `Self::default()`, all elements are equal to zero.
    /// 
    /// ```
    /// # use scilib::coordinate::Cartesian;
    /// let m = Cartesian { x: 0.0, y: 0.0, z: 0.0 };
    /// let n = Cartesian::new();
    /// let d = Cartesian::default();
    /// 
    /// assert_eq!(m, n);
    /// assert_eq!(n, d);
    /// ```
    pub fn new() -> Self {
        Self::default()
    }

    /// # From three points
    /// 
    /// Creates a Cartesian struct from three given points in space.
    /// 
    /// ```
    /// # use scilib::coordinate::Cartesian;
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

    /// # Distance between two points
    /// 
    /// Computes the distance between two points in space.
    pub fn distance(&self, other: Self) -> f64 {

        ((other.x - self.x).powi(2) + (other.y - self.y).powi(2) + (other.z - self.z).powi(2)).sqrt()
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # Addition
/// 
/// ```
/// # use scilib::coordinate::Cartesian;
/// let c1 = Cartesian::from(1.0, 2.5, 3.0);
/// let c2 = Cartesian::from(-1.0, 1, 0.0);
/// let res = c1 + c2;
/// let expected = Cartesian::from(0.0, 3.5, 3.0);
/// 
/// assert_eq!(res, expected);
/// ```
impl Add for Cartesian {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
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
/// # use scilib::coordinate::Cartesian;
/// let mut c1 = Cartesian::from(1.0, 2.5, 3.0);
/// let c2 = Cartesian::from(-1.0, 1, 0.0);
/// c1 += c2;
/// let expected = Cartesian::from(0.0, 3.5, 3.0);
/// 
/// assert_eq!(c1, expected);
/// ```
impl AddAssign for Cartesian {
    fn add_assign(&mut self, rhs: Self) {
        self.x += rhs.x;
        self.y += rhs.y;
        self.z += rhs.z;
    }
}

/// # Subtraction
/// 
/// ```
/// # use scilib::coordinate::Cartesian;
/// let c1 = Cartesian::from(1.0, 2.5, 3.0);
/// let c2 = Cartesian::from(-1.0, 1, 0.0);
/// let res = c1 - c2;
/// let expected = Cartesian::from(2, 1.5, 3.0);
/// 
/// assert_eq!(res, expected);
/// ```
impl Sub for Cartesian {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
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
/// # use scilib::coordinate::Cartesian;
/// let mut c1 = Cartesian::from(1.0, 2.5, 3.0);
/// let c2 = Cartesian::from(-1.0, 1, 0.0);
/// c1 -= c2;
/// let expected = Cartesian::from(2, 1.5, 3.0);
/// 
/// assert_eq!(c1, expected);
impl SubAssign for Cartesian {
    fn sub_assign(&mut self, rhs: Self) {
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
/// # use scilib::coordinate::Cartesian;
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
/// # use scilib::coordinate::Cartesian;
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
/// # use scilib::coordinate::Cartesian;
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
/// # use scilib::coordinate::Cartesian;
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
/// # use scilib::coordinate::Cartesian;
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
