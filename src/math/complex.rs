//!
//! # Complex numbers
//! 
//! The complex numbers module is secondary to other objectives of the crate. The functionalities will be
//! added to match the need of the other function, such as spherical harmonics, or Bessel functions.
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
    DivAssign       // Assigning division
};

use std::fmt::{     // Formatter display
    Display,        // The display itself
    Result as DRes  // The associated result
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # Complex structure
/// 
/// The principle is simple, we create both parts in a struct and treat them accordingly.
#[derive(Debug, Default)]
pub struct Complex {
    /// The real part of the number
    pub re: f64,
    /// The imaginary part of the number
    pub im: f64
}

/// # Display
/// 
/// Simply display the value of the number as a String
impl Display for Complex {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> DRes {
        write!(f, "{}", format!("{} {:+}i", self.re, self.im))?;
        Ok(())
    }
}

/// Implementing required methods
impl Complex {
    /// # From any numbers
    /// 
    /// Both parts can be any number that can be cast to `f64`.
    /// ```
    /// # use abacus::math::complex::Complex;
    /// let c1 = Complex::from(3.0, 2.0);
    /// let c2 = Complex::from(10, 1.5);
    /// let c3 = Complex::from(1, -2);
    /// 
    /// // All numbers are now f64
    /// assert!(c1.re == 3.0 && c1.im == 2.0);
    /// assert!(c2.re == 10.0 && c2.im == 1.5);
    /// assert!(c3.re == 1.0 && c3.im == -2.0);
    /// ```
    pub fn from<T, U>(re: T, im: U) -> Self
    where T: Into<f64>, U: Into<f64> {
        Self {
            re: re.into(),
            im: im.into()
        }
    }

    /// # Pure complex unity
    /// 
    /// Simply returns i.
    /// 
    /// ```
    /// # use abacus::math::complex::Complex;
    /// let res = Complex::i();
    /// 
    /// assert!(res.re == 0.0 && res.im ==1.0);
    /// ```
    pub fn i() -> Self {
        Self {
            re: 0.0,
            im: 1.0
        }
    }

    /// # Modulus computation
    /// 
    /// The modulus of a complex number is defined as the square root of the
    /// sum of its squared part.
    /// 
    /// ```
    /// # use abacus::math::complex::Complex;
    /// let c1 = Complex::from(2, -1.5);
    /// let c2 = Complex::from(-5.1, 17);
    /// 
    /// assert_eq!(c1.modulus(), 2.5);
    /// assert!((c2.modulus() - 17.7485210651).abs() < 1.0e-8);
    /// ```
    pub fn modulus(&self) -> f64 {
        (self.re.powi(2) + self.im.powi(2)).sqrt()
    }

    /// # Exponential
    /// 
    /// Computes the exponential value of a complex number
    /// 
    /// ```
    /// # use abacus::math::complex::Complex;
    /// let c = Complex::from(2, 2.2);
    /// let res = c.exp();
    /// 
    /// assert!((res.re - -4.3484677696).abs() < 1.0e-10);
    /// assert!((res.im - 5.97402528360).abs() < 1.0e-10);
    /// ```
    pub fn exp(&self) -> Self {
        let e: f64 = self.re.exp();
        Self {
            re: self.im.cos() * e,
            im: self.im.sin() * e
        }
    }

    /// # Complex conjugation
    /// 
    /// Conjugating a complex number changes the sign of the imaginary part.
    /// 
    /// ```
    /// # use abacus::math::complex::Complex;
    /// let c = Complex::from(3, 4.6);
    /// let c_star = c.conjugate();
    /// 
    /// assert!(c_star.re == 3.0 && c_star.im == -4.6);
    /// ```
    pub fn conjugate(&self) -> Self {
        Self {
            re: self.re,
            im: -self.im
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Here comes a long list of implementations for the operations

/// # Addition of complex numbers
/// 
/// ```
/// # use abacus::math::complex::Complex;
/// let c1 = Complex::from(2.1, 3.0);
/// let c2 = Complex::from(5.0, 0.5);
/// let res = c1 + c2;
/// 
/// assert!(res.re == 7.1 && res.im == 3.5);
/// ```
impl Add for Complex {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        Self {
            re: self.re + rhs.re,
            im: self.im + rhs.im
        }
    }
}

/// # Assigning Addition
/// 
/// ```
/// # use abacus::math::complex::Complex;
/// let mut c1 = Complex::from(2.1, 3.0);
/// let c2 = Complex::from(5.0, 0.5);
/// c1 += c2;
/// 
/// assert!(c1.re == 7.1 && c1.im == 3.5);
/// ```
impl AddAssign for Complex {
    fn add_assign(&mut self, rhs: Self) {
        self.re += rhs.re;
        self.im += rhs.im;
    }
}

/// # Subtraction
/// 
/// ```
/// # use abacus::math::complex::Complex;
/// let c1 = Complex::from(2.1, 3.0);
/// let c2 = Complex::from(5.0, 0.5);
/// let res = c1 - c2;
/// 
/// assert!(res.re == -2.9 && res.im == 2.5);
impl Sub for Complex {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        Self {
            re: self.re - rhs.re,
            im: self.im - rhs.im
        }
    }
}

/// # Assigning subtraction
/// 
/// ```
/// # use abacus::math::complex::Complex;
/// let mut c1 = Complex::from(2.1, 3.0);
/// let c2 = Complex::from(5.0, 0.5);
/// c1 -= c2;
/// 
/// assert!(c1.re == -2.9 && c1.im == 2.5);
impl SubAssign for Complex {
    fn sub_assign(&mut self, rhs: Self) {
        self.re -= rhs.re;
        self.im -= rhs.im;
    }
}

/// # Multiplication
/// 
/// ```
/// # use abacus::math::complex::Complex;
/// let c1 = Complex::from(2.1, 3.0);
/// let c2 = Complex::from(5.0, 0.5);
/// let res = c1 * c2;
/// 
/// assert!(res.re == 9.0 && res.im == 16.05);
impl Mul for Complex {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        Self {
            re: self.re * rhs.re - self.im * rhs.im,
            im: self.re * rhs.im + self.im * rhs.re
        }
    }
}

/// # Multiplication to f64 (real): `c * f64`
/// 
/// ```
/// # use abacus::math::complex::Complex;
/// let c = Complex::from(5, 2.0);
/// let res = c * 3.0;
/// 
/// assert!(res.re == 15.0 && res.im == 6.0);
impl Mul<f64> for Complex {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self::Output {
        Self {
            re: self.re * rhs,
            im: self.im * rhs
        }
    }
}

/// # Multiplication to f64 (real): `f64 * c`
/// 
/// ```
/// # use abacus::math::complex::Complex;
/// let c = Complex::from(5, 2.0);
/// let res = 3.0 * c;
/// 
/// assert!(res.re == 15.0 && res.im == 6.0);
impl Mul<Complex> for f64 {
    type Output = Complex;
    fn mul(self, rhs: Complex) -> Self::Output {
        Complex {
            re: self * rhs.re,
            im: self * rhs.im
        }
    }
}

/// # Assigning multiplication
/// 
/// ```
/// # use abacus::math::complex::Complex;
/// let mut c1 = Complex::from(2.1, 3.0);
/// let c2 = Complex::from(5.0, 0.5);
/// c1 *= c2;
/// 
/// assert!(c1.re == 9.0 && c1.im == 16.05);
impl MulAssign for Complex {
    fn mul_assign(&mut self, rhs: Self) {
        let old_re: f64 = self.re;
        self.re = self.re * rhs.re - self.im * rhs.im;
        self.im = old_re * rhs.im + self.im * rhs.re;
    }
}

/// # Division
/// 
/// ```
/// # use abacus::math::complex::Complex;
/// let c1 = Complex::from(2.1, 3.0);
/// let c2 = Complex::from(5.0, 0.5);
/// let res = c1 / c2;
/// 
/// assert!((res.re - 0.47524752475).abs() < 1.0e-9 && (res.im - 0.5524752475).abs() < 1.0e-9);
impl Div for Complex {
    type Output = Self;
    fn div(self, rhs: Self) -> Self::Output {
        let div: f64 = rhs.re.powi(2) + rhs.im.powi(2);
        Self {
            re: (self.re * rhs.re + self.im * rhs.im) / div,
            im: (self.im * rhs.re - self.re * rhs.im) / div
        }
    }
}

/// # Assigning division
/// 
/// ```
/// # use abacus::math::complex::Complex;
/// let mut c1 = Complex::from(2.1, 3.0);
/// let c2 = Complex::from(5.0, 0.5);
/// c1 /= c2;
/// 
/// assert!((c1.re - 0.47524752475).abs() < 1.0e-9 && (c1.im - 0.5524752475).abs() < 1.0e-9);
impl DivAssign for Complex {
    fn div_assign(&mut self, rhs: Self) {
        let div: f64 = rhs.re.powi(2) + rhs.im.powi(2);
        let old_re: f64 = self.re;
        self.re = (self.re * rhs.re + self.im * rhs.im) / div;
        self.im = (self.im * rhs.re - old_re * rhs.im) / div;
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
