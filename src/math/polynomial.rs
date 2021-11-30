//!
//! Classical polynomials
//! 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

use super::{            // Using parts from the crate
    basic,              // Basic functions
    complex::Complex    // Using Complex numbers
};

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

/// # Bernoulli polynomials
#[derive(Debug, Default)]
pub struct Bernoulli {
    /// The order of the polynomial
    pub n: usize,
    /// The factor of each term
    factor: Vec<f64>,
    /// The power at each term
    power: Vec<i32>
}

/// Display for the Laguerre polynomials
impl std::fmt::Display for Bernoulli {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        
        let mut s: String = String::from("");

        for (f, p) in self.factor.iter().zip(&self.power) {
            match p {
                0 => s += &format!("{:+} ", f),
                1 => s += &format!("{:+}x ", f),
                _ => s += &format!("{:+}x^{} ", f, p)
            }
            
        }
        write!(f, "[n={}] :: {}", self.n, s)?;
        Ok(())
    }
}

/// Implementing the required methods
impl Bernoulli {

    /// Bernoulli number generator
    /// 
    /// Following the Bm- convention.
    /// ```
    /// # use scilib::math::polynomial::Bernoulli;
    /// let num_0: f64 = Bernoulli::gen_number(0);
    /// let num_1: f64 = Bernoulli::gen_number(1);
    /// let num_2: f64 = Bernoulli::gen_number(2);
    /// 
    /// assert!((num_0 - 1.0).abs() <= 1.0e-8);
    /// assert!((num_1 - -0.5).abs() <= 1.0e-8);
    /// assert!((num_2 - 1.0 / 6.0).abs() <= 1.0e-8);
    /// ```
    pub fn gen_number(m: usize) -> f64 {

        let mut res: f64 = 0.0;     // Resulting number
        
        // First loop
        for k in 0..=m {
            // Precomputing the factor
            let k1: f64 = (k + 1) as f64;
            // Inner loop
            for v in 0..=k {
                let sign: f64 = (-1.0_f64).powi(v as i32);      // -1^v
                let binom: f64 = basic::binomial(k, v) as f64;  // The binomial value
                let frac: f64 = (v as f64).powi(m as i32) / k1; // The fraction

                res += sign * binom * frac;                     // Summing the current value
            }
        }

        // Returning the result
        res
    }

    /// Produces the factors and powers for nth order polynomial.
    /// 
    /// Returns: `Self`, the corresponding struct
    /// ```
    /// # use scilib::math::polynomial::Bernoulli;
    /// let p2 = Bernoulli::new(2);         // n=2
    /// let p3 = Bernoulli::new(3);         // n=3
    /// ```
    pub fn new(n: usize) -> Self {

        // Initializing the vectors
        let mut factor: Vec<f64> = Vec::new();
        let mut power: Vec<i32> = Vec::new();

        for k in 0..=n {
            power.push(k as i32);
            let coef: f64 = basic::binomial(n, k) as f64 * Self::gen_number(n - k);
            factor.push(coef);
        }

        // Returning associated struct
        Self {
            n,
            factor,
            power
        }
    }

    /// Computes the value of `x` for the given polynomial.
    /// 
    /// Returns: the result of the polynomial Bn(x)
    /// 
    /// ```
    /// # use scilib::math::polynomial::Bernoulli;
    /// let x = 2.5;                        // Example value
    /// 
    /// let p = Bernoulli::new(3);          // n=3
    /// 
    /// // Computing the results for the polynomial
    /// let res = p.compute(x);
    ///
    /// // Comparing to tabulated values
    /// assert!((res - 7.5).abs() <= 1.0e-8);
    /// ```
    pub fn compute(&self, x: f64) -> f64 {
        // Iterates through the values of the factors and powers
        self.factor.iter().zip(&self.power).fold(0.0, |res, (f, p)| res + f * x.powi(*p))
    }

    /// Computes the value of `z` for the given polynomial (z: complex).
    /// 
    /// Returns: the result of the polynomial Bn(z)
    pub fn compute_complex(&self, z: Complex) -> Complex {
        // Iterates through the values of the factors and powers
        self.factor.iter().zip(&self.power).fold(Complex::unity(), |res, (f, p)| res + *f * z.powi(*p))
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # Bernoulli polynomials
#[derive(Debug, Default)]
pub struct Euler {
    /// The order of the polynomial
    pub n: usize,
    /// The factor of each term
    factor: Vec<f64>,
    /// The power at each term
    power: Vec<i32>
}

/// Display for the Laguerre polynomials
impl std::fmt::Display for Euler {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        
        let mut s: String = String::from("");

        for (f, p) in self.factor.iter().zip(&self.power) {
            match p {
                0 => s += &format!("{:+} ", f),
                1 => s += &format!("{:+}x ", f),
                _ => s += &format!("{:+}x^{} ", f, p)
            }
            
        }
        write!(f, "[n={}] :: {}", self.n, s)?;
        Ok(())
    }
}

impl Euler {
    /// Euler number generator
    /// 
    /// ```
    /// # use scilib::math::polynomial::Euler;
    /// // We know that all odd Euler numbers are 0
    /// for odd in vec![1, 3, 5, 7, 9, 11] {
    ///     let ek: f64 = Euler::gen_number(odd);
    ///     assert_eq!(ek, 0.0);
    /// }
    /// 
    /// // The others oscillate between signs
    /// let e0: f64 = Euler::gen_number(0);
    /// let e2: f64 = Euler::gen_number(2);
    /// let e10: f64 = Euler::gen_number(10);
    /// 
    /// assert_eq!(e0, 1.0);
    /// assert_eq!(e2, -1.0);
    /// assert_eq!(e10, -50_521.0);
    /// ```
    pub fn gen_number(m: usize) -> f64 {

        let mut res: f64 = 0.0;
        
        // All odd Euler numbers are 0
        if m % 2 != 0 {
            return res;
        } else if m == 0 {
            return 1.0;
        }

        for i in 1..=m {
            let sign: f64 = (-1.0_f64).powi(i as i32);
            let fact: f64 = 1.0 / 2.0_f64.powi(i as i32);
            let mut term: f64 = 0.0;

            for l in 0..=(2*i) {
                let sign2: f64 = (-1.0_f64).powi(l as i32);
                let binom: f64 = basic::binomial(2 * i, l) as f64;
                let par: f64 = (i as f64 - l as f64).powi(m as i32);

                term += sign2 * binom * par;
            }

            res += sign * fact * term;
        }

        res
    }

    /// Produces the factors and powers for nth order polynomial.
    /// 
    /// Returns: `Self`, the corresponding struct
    /// ```
    /// # use scilib::math::polynomial::Euler;
    /// let p3 = Euler::new(3);             // n=3
    /// let p7 = Euler::new(7);             // n=7
    /// ```
    pub fn new(n: usize) -> Self {

        // Initializing the vectors
        let mut factor: Vec<f64> = vec![0.0; n + 1];
        let power: Vec<i32> = (0..=n).map(|x| x as i32).collect();

        for k in 0..=n {
            let binom: f64 = basic::binomial(n, k) as f64;
            let f: f64 = Self::gen_number(k) / 2.0_f64.powi(k as i32);

            // Second loop from the (x - 1/2)^(n-k)
            for p in 0..=(n-k) {
                let pre: f64 = (-0.5_f64).powi(p as i32);
                let triangle_val: f64 = basic::binomial(n-k, p) as f64;

                factor[n-k-p] += binom * f * pre * triangle_val;
            }
        }

        // Returning associated struct
        Self {
            n,
            factor,
            power
        }
    }

    /// Computes the value of `x` for the given polynomial (x: real).
    /// 
    /// Returns: the result of the polynomial En(x)
    /// 
    /// ```
    /// # use scilib::math::polynomial::Euler;
    /// let x = -1.1;                    // Example value
    /// 
    /// let p = Euler::new(5);          // n=5
    /// 
    /// // Computing the results for the polynomial
    /// let res = p.compute(x);
    ///
    /// // Comparing to tabulated values
    /// assert!((res - -2.74576).abs() <= 1.0e-8);
    /// ```
    pub fn compute(&self, x: f64) -> f64 {
        // Iterates through the values of the factors and powers
        self.factor.iter().zip(&self.power).fold(0.0, |res, (f, p)| res + f * x.powi(*p))
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
