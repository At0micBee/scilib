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
    DivAssign,      // Assigning division
    Neg             // Negation
};

use std::fmt::{     // Formatter display
    Display,        // The display itself
    Result as DRes  // The associated result
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # Complex structure
/// 
/// The principle is simple, we create both parts in a struct and treat them accordingly.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Complex {
    /// The real part of the number
    pub re: f64,
    /// The imaginary part of the number
    pub im: f64
}

/// # Display
/// 
/// Returns the complex in the for a + bi, whe the sign of b is always showing.
impl Display for Complex {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> DRes {
        write!(f, "{}", format!("{} {:+}i", self.re, self.im))?;
        Ok(())
    }
}

/// Implementing required methods
impl Complex {
    /// # New Complex
    /// 
    /// Simply returns 0 +0i.
    /// 
    /// ```
    /// # use scilib::math::complex::Complex;
    /// let res = Complex::new();
    /// 
    /// assert!(res.re == 0.0 && res.im == 0.0);
    /// ```
    pub const fn new() -> Self {
        Self {
            re: 0.0,
            im: 0.0
        }
    }

    /// # Pure complex unity
    /// 
    /// Simply returns 0 +1i.
    /// 
    /// ```
    /// # use scilib::math::complex::Complex;
    /// let res = Complex::i();
    /// 
    /// assert!(res.re == 0.0 && res.im == 1.0);
    /// ```
    pub const fn i() -> Self {
        Self {
            re: 0.0,
            im: 1.0
        }
    }

    /// # Pure real unity
    /// 
    /// Simply returns 1 + 0i.
    /// 
    /// ```
    /// # use scilib::math::complex::Complex;
    /// let res = Complex::unity();
    /// 
    /// assert!(res.re == 1.0 && res.im == 0.0);
    /// ```
    pub const fn unity() -> Self {
        Self {
            re: 1.0,
            im: 0.0
        }
    }

    /// # From any numbers
    /// 
    /// Both parts can be any number that can be cast to `f64`.
    /// ```
    /// # use scilib::math::complex::Complex;
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

    /// # From polar coordinates
    /// 
    /// Creates the complex number based on polar coordinates values
    /// 
    /// ```
    /// # use std::f64::consts::PI;
    /// # use std::f64::consts::SQRT_2;
    /// # use scilib::math::complex::Complex;
    /// let c = Complex::from_polar(PI / 4.0, SQRT_2);
    /// 
    /// assert!((c.re - 1.0).abs() < 1.0e-8 && c.im == 1.0);
    /// ```
    pub fn from_polar<T, U>(arg: T, norm: U) -> Self
    where T: Into<f64> + Copy, U: Into<f64> + Copy {
        Self::from(arg.into().cos() * norm.into(), arg.into().sin() * norm.into())
    }

    /// # Exponential
    /// 
    /// Computes the exponential value of a complex number
    /// 
    /// ```
    /// # use scilib::math::complex::Complex;
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

    /// # Natural logarithm
    /// 
    /// Computes the `ln` of self.
    /// 
    /// ```
    /// # use scilib::math::complex::Complex;
    /// let c = Complex::from(1.2, 5.35);
    /// let res = c.ln();
    /// 
    /// assert!((res.re - 1.70163927393298).abs() < 1.0e-10);
    /// assert!((res.im - 1.35014914413455).abs() < 1.0e-10);
    /// ```
    pub fn ln(&self) -> Self {
        let (arg, norm): (f64, f64) = self.polar();
        Self {
            re: norm.ln(),
            im: arg
        }
    }

    /// # Cosine function
    /// 
    /// Computes the cosine value of the given complex number.
    /// Formula: `cos(a + ib) = cos(a)cosh(b) - i sin(a)sinh(b)`.
    /// 
    /// ```
    /// # use scilib::math::complex::Complex;
    /// let c = Complex::from(3, 2.1);
    /// let res = c.cos();
    /// 
    /// assert!((res.re - -4.102838942268).abs() < 1.0e-10);
    /// assert!((res.im - -0.567564455868).abs() < 1.0e-10);
    /// ```
    pub fn cos(&self) -> Self {
        Self {
            re: self.re.cos() * self.im.cosh(),
            im: -self.re.sin() * self.im.sinh()
        }
    }

    /// # Sinus function
    /// 
    /// Computes the sinus value of the given complex number.
    /// Formula: `sin(a + ib) = sin(a)cosh(b) + i cos(a)sinh(b)`
    /// 
    /// ```
    /// # use scilib::math::complex::Complex;
    /// let c = Complex::from(3, 2.1);
    /// let res = c.sin();
    /// 
    /// assert!((res.re - 0.5848455080109).abs() < 1.0e-10);
    /// assert!((res.im - -3.9816079971376).abs() < 1.0e-10);
    /// ```
    pub fn sin(&self) -> Self {
        Self {
            re: self.re.sin() * self.im.cosh(),
            im: self.re.cos() * self.im.sinh()
        }
    }

    /// # Tangent function
    /// 
    /// Computes the tangent value of the given complex number.
    /// Formula: `tan(x) = sin(x) / cos(x)`.
    /// 
    /// ```
    /// # use scilib::math::complex::Complex;
    /// let c = Complex::from(3, 2.1);
    /// let res = c.tan();
    /// 
    /// assert!((res.re - -0.0081436522788).abs() < 1.0e-10);
    /// assert!((res.im - 0.97157848523977).abs() < 1.0e-10);
    /// ```
    pub fn tan(&self) -> Self {
        self.sin() / self.cos()
    }

    /// # Complex conjugation
    /// 
    /// Conjugating a complex number changes the sign of the imaginary part.
    /// 
    /// ```
    /// # use scilib::math::complex::Complex;
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

    /// # Argument for polar coordinates
    /// 
    /// Computes the argument in the polar plan.
    /// 
    /// ```
    /// # use std::f64::consts::PI;
    /// # use scilib::math::complex::Complex;
    /// let c1 = Complex::from(1.0, 1.0);
    /// let c2 = Complex::from(2.5, -1);
    /// 
    /// assert!((c1.arg() - PI/4.0).abs() < 1.0e-8);
    /// assert!((c2.arg() - -0.3805063771).abs() < 1.0e-8);
    /// ```
    pub fn arg(&self) -> f64 {
        self.im.atan2(self.re)
    }

    /// # Modulus computation
    /// 
    /// The modulus of a complex number is defined as the square root of the
    /// sum of its squared part.
    /// 
    /// ```
    /// # use scilib::math::complex::Complex;
    /// let c1 = Complex::from(2, -1.5);
    /// let c2 = Complex::from(-5.1, 17);
    /// 
    /// assert_eq!(c1.modulus(), 2.5);
    /// assert!((c2.modulus() - 17.7485210651).abs() < 1.0e-8);
    /// ```
    pub fn modulus(&self) -> f64 {
        (self.re.powi(2) + self.im.powi(2)).sqrt()
    }

    /// # The polar coordinates of the number
    /// 
    /// Returns a tuple where the zeroth element is the argument and the first
    /// element is the modulus (or norm) of the number.
    /// 
    /// ```
    /// # use std::f64::consts::PI;
    /// # use std::f64::consts::SQRT_2;
    /// # use scilib::math::complex::Complex;
    /// let c = Complex::from(1, 1);
    /// let p = c.polar();
    /// 
    /// assert!((p.0 - PI/4.0).abs() < 1.0e-8);
    /// assert!((p.1 - SQRT_2).abs() < 1.0e-8);
    /// ```
    pub fn polar(&self) -> (f64, f64) {
        (self.arg(), self.modulus())
    }

    /// # Raising to an integer power
    ///
    /// ```
    /// # use scilib::math::complex::Complex;
    /// let c = Complex::from(3.0, 1.0);
    /// let res1 = c.powi(4);
    /// let res2 = c.powi(-1);
    /// 
    /// assert!(res1.re == 28.0 && res1.im == 96.0);
    /// assert!(res2.re == 0.3 && res2.im == -0.1);
    /// ```
    pub fn powi(self, e: i32) -> Self {

        // x^0 = 1, technically not correct for 0^0 but we have to do something
        if e == 0 {
            return 1.0.into();
        }

        let mut res: Self = self;

        // When e > 0
        for _ in 1..e {
            res *= self;
        }

        // When e < 0
        for _ in e..=0 {
            res /= self;
        }

        res
    }

    /// # Raising to a real power
    ///
    /// ```
    /// # use scilib::math::complex::Complex;
    /// let c = Complex::from(2.5, -3.6);
    /// let res = c.powf(1.25);
    /// 
    /// assert!((res.re - 2.2697926495).abs() < 1.0e-8 && (res.im - -5.9215705908).abs() < 1.0e-8);
    /// ```
    pub fn powf(&self, e: f64) -> Self {
        // Cheating via polar
        let (arg, norm): (f64, f64) = self.polar();
        Self::from_polar(arg * e, norm.powf(e))
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Here comes a long list of implementations for the operations

/// # Conversion from a scalar
/// 
/// Takes a scalar value and assigns it to the real part, as long as the type
/// allows conversion to `f64`.
/// 
/// ```
/// # use scilib::math::complex::Complex;
/// let c1: Complex = 3.5.into();
/// let c2: Complex = (-5).into();
/// 
/// assert!(c1.re == 3.5 && c1.im == 0.0);
/// assert!(c2.re == -5.0 && c2.im == 0.0);
/// ```
impl<T: Into<f64>> From<T> for Complex {
    fn from(val: T) -> Self {
        Self {
            re: val.into(),
            im: 0.0
        }
    }
}

/// # Addition of complex numbers
/// 
/// ```
/// # use scilib::math::complex::Complex;
/// let c1 = Complex::from(2.1, 3.0);
/// let c2 = Complex::from(5.0, 0.5);
/// let res = c1 + c2;
/// let res2 = c1 + 5.5;
/// 
/// assert!(res.re == 7.1 && res.im == 3.5);
/// assert!(res2.re == 7.6 && res2.im == 3.0);
/// ```
impl<T: Into<Self>> Add<T> for Complex {
    type Output = Self;
    fn add(self, rhs: T) -> Self::Output {
        let rhs: Self = rhs.into();
        Self {
            re: self.re + rhs.re,
            im: self.im + rhs.im
        }
    }
}

/// # Addition to f64 (real): `f64 + c`
/// 
/// ```
/// # use scilib::math::complex::Complex;
/// let c = Complex::from(7, 2.0);
/// let res = 3.0 + c;
/// 
/// assert!(res.re == 10.0 && res.im == 2.0);
impl Add<Complex> for f64 {
    type Output = Complex;
    fn add(self, rhs: Complex) -> Self::Output {
        Complex {
            re: self + rhs.re,
            im: self + rhs.im
        }
    }
}

/// # Assigning Addition
/// 
/// ```
/// # use scilib::math::complex::Complex;
/// let mut c1 = Complex::from(2.1, 3.0);
/// let mut c2 = Complex::from(5.0, 0.5);
/// c1 += c2;
/// c2 += 2;
/// 
/// assert!(c1.re == 7.1 && c1.im == 3.5);
/// assert!(c2.re == 7.0 && c2.im == 0.5);
/// ```
impl<T: Into<Self>> AddAssign<T> for Complex {
    fn add_assign(&mut self, rhs: T) {
        let rhs: Self = rhs.into();
        self.re += rhs.re;
        self.im += rhs.im;
    }
}

/// # Subtraction
/// 
/// ```
/// # use scilib::math::complex::Complex;
/// let c1 = Complex::from(2.1, 3.0);
/// let c2 = Complex::from(5.0, 0.5);
/// let res = c1 - c2;
/// let res2 = c1 - 5;
/// 
/// assert!(res.re == -2.9 && res.im == 2.5);
/// assert!(res2.re == -2.9 && res2.im == 3.0);
impl<T: Into<Self>> Sub<T> for Complex {
    type Output = Self;
    fn sub(self, rhs: T) -> Self::Output {
        let rhs: Self = rhs.into();
        Self {
            re: self.re - rhs.re,
            im: self.im - rhs.im
        }
    }
}

/// # Subtraction to f64 (real): `f64 - c`
/// 
/// ```
/// # use scilib::math::complex::Complex;
/// let c = Complex::from(10, 2.0);
/// let res = 3.0 - c;
/// 
/// assert!(res.re == 7.0 && res.im == 2.0);
impl Sub<Complex> for f64 {
    type Output = Complex;
    fn sub(self, rhs: Complex) -> Self::Output {
        Complex {
            re: self - rhs.re,
            im: self - rhs.im
        }
    }
}

/// # Assigning subtraction
/// 
/// ```
/// # use scilib::math::complex::Complex;
/// let mut c1 = Complex::from(2.1, 3.0);
/// let mut c2 = Complex::from(5.0, 0.5);
/// c1 -= c2;
/// c2 -= 12.0;
/// 
/// assert!(c1.re == -2.9 && c1.im == 2.5);
/// assert!(c2.re == -7.0 && c2.im == 0.5)
impl<T: Into<Self>> SubAssign<T> for Complex {
    fn sub_assign(&mut self, rhs: T) {
        let rhs: Self = rhs.into();
        self.re -= rhs.re;
        self.im -= rhs.im;
    }
}

/// # Multiplication
/// 
/// ```
/// # use scilib::math::complex::Complex;
/// let c1 = Complex::from(2.1, 3.0);
/// let c2 = Complex::from(5.0, 0.5);
/// let res = c1 * c2;
/// let res2 = c1 * 2.0;
/// 
/// assert!(res.re == 9.0 && res.im == 16.05);
/// assert!(res2.re == 4.2 && res2.im == 6.0);
impl<T: Into<Self>> Mul<T> for Complex {
    type Output = Self;
    fn mul(self, rhs: T) -> Self::Output {
        let rhs: Self = rhs.into();
        Self {
            re: self.re * rhs.re - self.im * rhs.im,
            im: self.re * rhs.im + self.im * rhs.re
        }
    }
}

/// # Multiplication to f64 (real): `f64 * c`
/// 
/// ```
/// # use scilib::math::complex::Complex;
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
/// # use scilib::math::complex::Complex;
/// let mut c1 = Complex::from(2.1, 3.0);
/// let mut c2 = Complex::from(5.0, 0.5);
/// c1 *= c2;
/// c2 *= 2;
/// 
/// assert!(c1.re == 9.0 && c1.im == 16.05);
/// assert!(c2.re == 10.0 && c2.im == 1.0);
impl<T: Into<Self>> MulAssign<T> for Complex {
    fn mul_assign(&mut self, rhs: T) {
        let rhs: Self = rhs.into();
        let old_re: f64 = self.re;
        self.re = self.re * rhs.re - self.im * rhs.im;
        self.im = old_re * rhs.im + self.im * rhs.re;
    }
}

/// # Division
/// 
/// ```
/// # use scilib::math::complex::Complex;
/// let c1 = Complex::from(2.1, 3.0);
/// let c2 = Complex::from(5.0, 0.5);
/// let res = c1 / c2;
/// let res2 = c1 / 2.0;
/// 
/// assert!((res.re - 0.47524752475).abs() < 1.0e-9 && (res.im - 0.5524752475).abs() < 1.0e-9);
/// assert!(res2.re == 1.05 && res2.im == 1.5);
impl<T: Into<Self>> Div<T> for Complex {
    type Output = Self;
    fn div(self, rhs: T) -> Self::Output {
        let rhs: Self = rhs.into();
        let div: f64 = rhs.re.powi(2) + rhs.im.powi(2);
        Self {
            re: (self.re * rhs.re + self.im * rhs.im) / div,
            im: (self.im * rhs.re - self.re * rhs.im) / div
        }
    }
}

/// # Division to f64 (real): `f64 / c`
/// 
/// ```
/// # use scilib::math::complex::Complex;
/// let c = Complex::from(2.0, 4.0);
/// let res = 2.0 / c;
/// 
/// assert!(res.re == 0.2 && res.im == 0.4);
impl Div<Complex> for f64 {
    type Output = Complex;
    fn div(self, rhs: Complex) -> Self::Output {
        // num / (a + ib) == a * num / (a² + b²) - i * b * num / (a² + b²)
        let modulus_squared = rhs.re * rhs.re + rhs.im * rhs.im;
        Complex {
            re: self * rhs.re / modulus_squared,
            im: -self * rhs.im / modulus_squared
        }
    }
}

/// # Assigning division
/// 
/// ```
/// # use scilib::math::complex::Complex;
/// let mut c1 = Complex::from(2.1, 3.0);
/// let mut c2 = Complex::from(5.0, 0.5);
/// c1 /= c2;
/// c2 /= 10.0;
/// 
/// assert!((c1.re - 0.47524752475).abs() < 1.0e-9 && (c1.im - 0.5524752475).abs() < 1.0e-9);
/// assert!(c2.re == 0.5 && c2.im == 0.05);
impl<T: Into<Self>> DivAssign<T> for Complex {
    fn div_assign(&mut self, rhs: T) {
        let rhs: Self = rhs.into();
        let div: f64 = rhs.re.powi(2) + rhs.im.powi(2);
        let old_re: f64 = self.re;
        self.re = (self.re * rhs.re + self.im * rhs.im) / div;
        self.im = (self.im * rhs.re - old_re * rhs.im) / div;
    }
}

/// # Negation
/// 
/// Returns the opposite of the number
/// 
/// ```
/// # use scilib::math::complex::Complex;
/// let c = Complex::from(1, 0.05);
/// let c_neg = -c;
/// 
/// assert!(c_neg.re == -1.0 && c_neg.im ==-0.05);
/// ```
impl Neg for Complex {
    type Output = Self;
    fn neg(self) -> Self::Output {
        Self {
            re: -self.re,
            im: -self.im
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
