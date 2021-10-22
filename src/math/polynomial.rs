//!
//! Classical polynomials
//! 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

use super::basic;   // Basic functions

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # Legendre polynomials
#[derive(Debug, Default)]
pub struct Legendre {
    /// The order of the polynomial
    pub l: usize,
    /// The derivative order
    pub m: i32,
    /// The factor of each term
    factor: Vec<f64>,
    /// The power at each term
    power: Vec<i32>,
    /// The pre-factor associated to m
    pre_f: f64
}

/// Display for the Legendre polynomials
impl std::fmt::Display for Legendre {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        
        let mut s: String = String::from(format!("{} * ( ", self.pre_f));

        for (f, p) in self.factor.iter().zip(&self.power) {
            match p {
                0 => s += &format!("{:+} ", f),
                1 => s += &format!("{:+}x ", f),
                _ => s += &format!("{:+}x^{} ", f, p)
            }
            
        }
        write!(f, "[l={}, m={}] :: {})", self.l, self.m, s)?;
        Ok(())
    }
}

/// Implementing the required methods
impl Legendre {
    /// Produces the factors and powers for nth order polynomial, where
    /// `l` is the order of the Legendre polynomial and `m` is the derivative order.
    /// 
    /// Returns: `Self`, the corresponding struct
    ///
    /// ```
    /// # use scilib::math::polynomial::Legendre;
    /// let p20 = Legendre::new(2, 0);      // l=2, m=0
    /// let p21 = Legendre::new(2, 1);      // l=2, m=1
    /// ```
    pub fn new(l: usize, m: i32) -> Self {

        // Checking that the range is good
        assert!(m >= -(l as i32) && m <= l as i32);

        // Initializing the vectors
        let mut factor: Vec<f64> = Vec::new();
        let mut power: Vec<i32> = Vec::new();
        
        // Going through the powers of the order
        for k in 0..=(l / 2) {
            power.push((l - 2 * k) as i32);
            let coef: f64 = ((-1.0_f64).powi(k as i32) * (basic::binomial(l, k) * basic::binomial(2*l - 2*k, l)) as f64) / 2.0_f64.powi(l as i32);
            factor.push(coef);
        }

        // Computing the pre-factor associated to m
        let mut pre_f: f64 = (-1.0_f64).powi(m);    // Condon–Shortley phase
        if m < 0 {
            pre_f *= (-1_f64).powi(m) * (basic::factorial((l as i32 - m) as usize) / basic::factorial((l as i32 + m) as usize)) as f64;
        }

        // Returning associated struct
        let mut poly: Self = Self {
            l,
            m,
            factor,
            power,
            pre_f
        };

        // We derive the polynomial m times
        poly.derive(m.abs());

        // Returning the final polynomial
        poly
    }

    /// Computes the m derivative of the polynomial
    ///
    /// This is used to produce the Pnm(x) version of the polynomial.
    fn derive(&mut self, m: i32) {

        // Creating the marker to know if we need to delete the last element
        // when the power reaches 0
        let mut marker: bool;
        
        // Looping the derivation, m times
        for _ in 1..=m {

            // Resetting at false, and starting the computation of the new terms
            marker = false;
            for (f, p) in self.factor.iter_mut().zip(&mut self.power) {
                match p {
                    0 => marker = true,     // If p is zero, we do not touch anything
                    _ => {                  // For anything else, we compute the derivative
                        *f *= *p as f64;    // Changing th factor value
                        *p -= 1;            // Changing the power value
                    }
                }
            }

            // If the marker is true, we remove the zero power, and the associated factor
            if marker {
                self.factor.pop();  // pop() method is fine, 0 can only exist once per derivative
                self.power.pop();   // and if so, it will always be at the last position
            }
        }
    }

    /// Computes the value of `x` for the given polynomial.
    /// 
    /// Returns: the result of the polynomial Plm(x)
    /// 
    /// ```
    /// # use scilib::math::polynomial::Legendre;
    /// let x = -0.25;                      // Example value
    /// 
    /// let p20 = Legendre::new(2, 0);      // l=2, m=0
    /// let p21 = Legendre::new(2, 1);      // l=2, m=1
    /// 
    /// // Computing the results for each polynomial
    /// let res20 = p20.compute(x);
    /// let res21 = p21.compute(x);
    ///
    /// // Comparing to tabulated values
    /// assert_eq!(res20, -0.40625);
    /// assert_eq!(res21, 0.7261843774138907);
    /// ```
    pub fn compute(&self, x: f64) -> f64 {
        // Iterates through the values of the factors and powers
        let pre: f64 = self.pre_f * (1.0 - x.powi(2)).powf(self.m as f64 / 2.0);
        pre * self.factor.iter().zip(&self.power).fold(0.0, |res, (f, p)| res + f * x.powi(*p))
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # Laguerre polynomials
#[derive(Debug, Default)]
pub struct Laguerre {
    /// The order of the polynomial
    pub l: usize,
    /// The derivative order
    pub m: i32,
    /// The factor of each term
    factor: Vec<f64>,
    /// The power at each term
    power: Vec<i32>
}

/// Display for the Laguerre polynomials
impl std::fmt::Display for Laguerre {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        
        let mut s: String = String::from("");

        for (f, p) in self.factor.iter().zip(&self.power) {
            match p {
                0 => s += &format!("{:+} ", f),
                1 => s += &format!("{:+}x ", f),
                _ => s += &format!("{:+}x^{} ", f, p)
            }
            
        }
        write!(f, "[l={}, m={}] :: {}", self.l, self.m, s)?;
        Ok(())
    }
}

/// Implementing the required methods
impl Laguerre {
    /// Produces the factors and powers for nth order polynomial, where
    /// `l`: the order of the Legendre polynomial and `m`: The derivative order.
    /// 
    /// Returns: `Self`, the corresponding struct
    /// ```
    /// # use scilib::math::polynomial::Laguerre;
    /// let p20 = Laguerre::new(2, 0);      // l=2, m=0
    /// let p21 = Laguerre::new(2, 1);      // l=2, m=1
    /// ```
    pub fn new(l: usize, m: i32) -> Self {

        // Checking that the range is good
        assert!(m >= 0 && m <= l as i32);

        // Initializing the vectors
        let mut factor: Vec<f64> = Vec::new();
        let mut power: Vec<i32> = Vec::new();

        // Going through the powers of the order
        for i in (0..=l).rev() {
            power.push(i as i32);
            let coef: f64 = (-1.0_f64).powi(i as i32) * basic::binomial(l + m as usize, l - i) as f64 / basic::factorial(i) as f64;
            factor.push(coef);
        }

        // Returning associated struct
        Self {
            l,
            m,
            factor,
            power
        }
    }

    /// Computes the value of `x` for the given polynomial.
    /// 
    /// Returns: the result of the polynomial Plm(x)
    /// 
    /// ```
    /// # use scilib::math::polynomial::Laguerre;
    /// let x = 0.2;                        // Example value
    /// 
    /// let p = Laguerre::new(2, 1);        // l=2, m=1
    /// 
    /// // Computing the results for the polynomial
    /// let res = p.compute(x);
    ///
    /// // Comparing to tabulated values
    /// assert_eq!(res, 2.42);
    /// ```
    pub fn compute(&self, x: f64) -> f64 {
        // Iterates through the values of the factors and powers
        self.factor.iter().zip(&self.power).fold(0.0, |res, (f, p)| res + f * x.powi(*p))
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
