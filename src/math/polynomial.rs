//!
//! Classical polynomials
//! 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

use std::collections::HashMap;  // Rust Hashmap

use super::basic;               // Basic functions

use num_complex::Complex64;     // Using complex numbers from the num crate

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # Polynomial implementation
#[derive(Clone)]
pub struct Poly {
    coef: HashMap<i32, f64>,
    l: Option<f64>,
    pre_f: Option<f64>,
    compute_fn: fn(&Self, f64) -> f64,
    compute_fnc: fn(&Self, Complex64) -> Complex64
}

/// Display for polynomials
impl std::fmt::Display for Poly {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        
        // We write the power and associated factors
        writeln!(f, "Power :: Factor")?;
        for (power, factor) in &self.coef {
            writeln!(f, "{:5} :: {}", power, factor)?;
        }

        // If there is a pre-factor we show it as well
        match self.pre_f {
            Some(pre) => writeln!(f, "Pre factor: {}", pre)?,
            None => ()
        }

        Ok(())
    }
}

// Default value for Poly
impl Default for Poly {
    fn default() -> Self {
        Self {
            coef: HashMap::new(),
            l: None,
            pre_f: None,
            compute_fn: Self::compute_base,
            compute_fnc: Self::compute_base_complex,
        }
    }
}

impl Poly {

    /// # Creates a new polynomial
    /// 
    /// ## Definition
    /// We create a polynomial with integer order powers, and real coefficients, such that:
    /// $$
    /// p = c_0 + c_1x + c_2x^2 + ... + c_nx^n
    /// $$
    /// 
    /// The computation of the result follows:
    /// $$
    /// r = c_ix^i
    /// $$
    /// where the $i$ indices are where non zero coefficients exist.
    /// 
    /// ## Inputs
    /// - `pow_fac`: the power and associated factor
    /// 
    /// ## Example
    /// ```
    /// # use scilib::math::polynomial::Poly;
    /// // 2 + x + 3x^2
    /// let p = Poly::from(&[(0, 2.0), (1, 1.0), (2, 3.0)]);
    /// let res = p.compute(2.0);
    /// assert_eq!(16.0, res);
    /// ```
    pub fn from(pow_fac: &[(i32, f64)]) -> Self {

        let coef: HashMap<i32, f64> = pow_fac.iter().map(|pf| *pf).collect();

        Self {
            coef,
            l: None,
            pre_f: None,
            compute_fn: Self::compute_base,
            compute_fnc: Self::compute_base_complex,
        }
    }

    //////////////////////////////////////////////////
    // Creating special polynomials

    /// # Legendre polynomials
    /// ## Definition
    /// The [Legendre polynomials](https://en.wikipedia.org/wiki/Legendre_polynomials) are used as solution to the
    /// Legendre differential equations, which can be written as:
    /// $$
    /// (1-x^2)\frac{d^2P_n(x)}{dx^2} - 2x\frac{dP_n(x)}{dx} + n(n+1)P_n(x) = 0
    /// $$
    /// 
    /// The close form used in the program is as follows:
    /// $$
    /// P_n(x) = \frac{1}{2^n}\sum_k^{\left\lfloor\frac{n}{2}\right\rfloor}(-1)^k\binom{n}{k}\binom{2n-2k}{n}x^{n-2k}
    /// $$
    /// 
    /// Which can be generalized with a given derivative order $l$, following the formula:
    /// $$
    /// P_n^l(x) = (-1)^l (1 - x^2)^{l/2} \frac{d^l P_n(x)}{dx^l}
    /// $$
    ///
    /// We can also compute the associated polynomials when $l<0$, following the convention:
    /// $$
    /// P_n^{-l}(x) = (-1)^l\frac{(n-l)!}{(n+l)!}P_n^l(x)
    /// $$
    /// 
    /// ## Inputs
    /// Produces the factors and powers for nth order polynomial, where
    /// - `n`: the order of the Legendre polynomial
    /// - `l`: the derivative order
    /// 
    /// By definition, we have that $n\ge0$ and $-n\le l\le n$.
    /// 
    /// Returns a `Poly` struct, corresponding to the associated Legendre polynomial.
    ///
    /// ## Example
    /// ```
    /// # use num_complex::Complex64;
    /// # use scilib::math::polynomial::Poly;
    /// let p71 = Poly::legendre(7, 1);     // n=7, l=1
    /// let p72 = Poly::legendre(7, -2);    // n=7, l=-2
    /// 
    /// let x = -0.25;                      // Example real
    /// let z = Complex64::new(-1.2, 0.2);  // Example complex
    /// 
    /// // Computing the results for each polynomial
    /// let res71 = p71.compute(x);
    /// let res72 = p72.compute_complex(z);
    /// 
    /// // Expected result
    /// let expected_c = Complex64::new(-0.1297952, -0.324460533333);
    ///
    /// // Comparing to tabulated values
    /// assert!((res71 - -0.681433146961).abs() < 1.0e-10);
    /// assert!((res72 - expected_c).norm() < 1.0e-10);
    /// ```
    pub fn legendre(n: usize, l: i32) -> Self {

        // Checking that the range is good
        assert!(l >= -(n as i32) && l <= n as i32, "The derivative order l isn't valid for the given n!");

        // Initializing the coefficients
        let mut coef: HashMap<i32, f64> = HashMap::new();
        
        // Going through the powers of the order
        for k in 0..=(n / 2) {
            let c: f64 = ((-1.0_f64).powi(k as i32) * (basic::binomial(n, k) * basic::binomial(2*n - 2*k, n)) as f64) / 2.0_f64.powi(n as i32);
            coef.insert((n - 2 * k) as i32, c);
        }

        // Computing the pre-factor associated to m
        let pre_f: f64 = if l >= 0 {
            (-1_f64).powi(l)
        } else {
            basic::factorial((n as i32 + l) as usize) as f64 / basic::factorial((n as i32 - l) as usize) as f64
        };

        // Returning associated struct
        let mut poly: Self = Self {
            coef,
            l: Some(l.abs() as f64),
            pre_f: Some(pre_f),
            compute_fn: Self::compute_legendre,
            compute_fnc: Self::compute_legendre_complex
        };

        // We derive the polynomial m times
        poly.derive(l.abs());

        // Returning the final polynomial
        poly
    }

    /// # Laguerre polynomials
    /// 
    /// ## Definition
    /// The [Laguerre polynomials](https://en.wikipedia.org/wiki/Laguerre_polynomials) are the solution to the Laguerre differential equation:
    /// $$
    /// x\frac{d^2L_n^{(\alpha)}}{dx^2} + (\alpha+1-x)\frac{dL_n^{(\alpha)}}{dx} + nL_n^{(\alpha)} = 0
    /// $$
    /// 
    /// In this crate, we use the generalized form of the Laguerre polynomial $L_n^{(\alpha)}$, using the closed form:
    /// $$
    /// L_n^{(\alpha)}(x) = \sum_{i=0}^{n} (-1)^i\binom{n+\alpha}{n-i}\frac{x^i}{i!}
    /// $$
    /// 
    /// Which yields the standard Laguerre polynomial for $\alpha=0$, but lets the polynomial be used
    /// to solve a wider variety of equations (such as radial wave function in quantum mechanics).
    /// The binomial is the generalized version, which allows real input. The factors of the
    /// polynomials are normalized.
    /// 
    /// ## Inputs
    /// Produces the factors and powers for nth order polynomial, where
    /// - `n`: the order of the Legendre polynomial ($n$)
    /// - `l`: the derivative order $\alpha$
    /// 
    /// Returns a `Poly` struct, corresponding to the associated Laguerre polynomial.
    /// 
    /// ## Example
    /// ```
    /// # use num_complex::Complex64;
    /// # use scilib::math::polynomial::Poly;
    /// let p21 = Poly::laguerre(2, 1.0);       // n=2, l=1
    /// let p73 = Poly::laguerre(7, 3.0);       // n=7, l=3
    /// let p515 = Poly::laguerre(5, -1.5);     // n=5, partial order l=-1.5
    /// 
    /// let x = 0.2;                            // Example real
    /// let z = Complex64::new(1.2, -0.4);      // Example complex
    /// 
    /// // Computing the results for the polynomial
    /// let res = p21.compute(x);
    /// let res_c = p73.compute_complex(z);
    /// let res_p = p515.compute_complex(z);
    /// 
    /// // Expected results
    /// let expected_c = Complex64::new(-7.429297330793, 10.152990394920);
    /// let expected_p = Complex64::new(0.310088583333, -0.058726333333);
    ///
    /// // Comparing to tabulated values
    /// assert_eq!(res, 2.42);
    /// assert!((res_c - expected_c).norm() < 1.0e-10);
    /// assert!((res_p - expected_p).norm() < 1.0e-10);
    /// ```
    pub fn laguerre<U>(n: usize, l: U) -> Self
    where U: Into<f64> {

        let alpha: f64 = l.into();

        // Initializing the vectors
        let mut coef: HashMap<i32, f64> = HashMap::new();

        // Going through the powers of the order
        for i in (0..=n).rev() {
            let c: f64 = (-1.0_f64).powi(i as i32) * basic::binomial_reduced(n as f64 + alpha, n - i) as f64 / basic::factorial(i) as f64;
            coef.insert(i as i32, c);
        }

        // Returning associated struct
        Self {
            coef,
            l: Some(alpha),
            pre_f: None,
            compute_fn: Self::compute_base,
            compute_fnc: Self::compute_base_complex
        }
    }

    /// # Bernoulli polynomials
    /// ## Definition
    /// The [Bernoulli polynomial](https://en.wikipedia.org/wiki/Bernoulli_polynomials) are present in a great variety
    /// of particular functions, and are closely related to the Euler polynomial (also in this crate).
    /// They are defined by the generating function:
    /// $$
    /// \frac{t\exp(xt)}{\exp(t) - 1} = \sum_{n=0}^{\infty}B_n(x)\frac{t^n}{n!}
    /// $$
    /// To generate the polynomials, we use the explicit closed form:
    /// $$
    /// B_n(x) = \sum_{k=0}^{n}\binom{n}{k}B_{n-k}x^{k}
    /// $$
    /// Where $B_{n-k}$ correspond to the $(n-k)^\mathrm{th}$ [Bernoulli number](https://en.wikipedia.org/wiki/Bernoulli_number),
    /// also available in this crate.
    /// 
    /// ## Inputs
    /// - `n`: the order of the polynomial
    /// 
    /// By definition, we have that $n\ge0$.
    /// 
    /// Returns a `Poly` struct, corresponding to the associated Laguerre polynomial.
    /// 
    /// ## Example
    /// ```
    /// # use num_complex::Complex64;
    /// # use scilib::math::polynomial::Poly;
    /// let p2 = Poly::bernoulli(2);        // n=2
    /// let p3 = Poly::bernoulli(3);        // n=3
    /// 
    /// let x = 2.5;                        // Example real
    /// let z = Complex64::new(-1.2, 0.2);  // Example complex
    /// 
    /// // Computing the results for the polynomial
    /// let res2 = p2.compute_complex(z);
    /// let res3 = p3.compute(x);
    /// 
    /// let expected_c = Complex64::new(2.766666666666, -0.68);
    ///
    /// // Comparing to tabulated values
    /// assert!((res3 - 7.5).abs() < 1.0e-8);
    /// assert!((res2 - expected_c).norm() < 1.0e-8);
    /// ```
    pub fn bernoulli(n: usize) -> Self {

        // Initializing the vectors
        let mut coef: HashMap<i32, f64> = HashMap::new();

        for k in 0..=n {
            let c: f64 = basic::binomial(n, k) as f64 * Self::bernoulli_number(n - k);
            coef.insert(k as i32, c);
        }

        // Returning associated struct
        Self {
            coef,
            l: None,
            pre_f: None,
            compute_fn: Self::compute_base,
            compute_fnc: Self::compute_base_complex
        }
    }

    /// # Euler polynomials
    /// ## Definition
    /// The [Euler polynomial](https://en.wikipedia.org/wiki/Bernoulli_polynomials) are present in a great variety
    /// of particular functions, and are closely related to the Euler polynomial (also in this crate).
    /// They are defined by the generating function:
    /// $$
    /// \frac{2\exp(xt)}{\exp(t) + 1} = \sum_{n=0}^{\infty}E_n(x)\frac{t^n}{n!}
    /// $$
    /// To generate the polynomials, we use the explicit formula:
    /// $$
    /// E_n(x) = \sum_{k=0}^{n}\binom{n}{k} \frac{E_k}{2^k}\left( x - \frac{1}{2} \right)^{n-k}
    /// $$
    /// Where $E_{k}$ correspond to the $k^\mathrm{th}$ [Euler number](https://en.wikipedia.org/wiki/Euler_numbers),
    /// also available in this crate.
    /// 
    /// ## Inputs
    /// - `n`: order of the polynomial
    /// 
    /// By definition, we have that $n\ge0$.
    /// 
    /// Returns a `Poly` struct, corresponding to the associated Laguerre polynomial.
    /// 
    /// ## Example
    /// ```
    /// # use num_complex::Complex64;
    /// # use scilib::math::polynomial::Poly;
    /// let p5 = Poly::euler(5);            // n=5
    /// let p6 = Poly::euler(6);            // n=7
    /// 
    /// let x = -1.1;                       // Example real
    /// let z = Complex64::new(1.0, -2.5);  // Example complex
    /// 
    /// // Computing the results for the polynomial
    /// let res5 = p5.compute(x);
    /// let res6 = p6.compute_complex(z);
    ///
    /// // Comparing to tabulated values
    /// assert!((res5 - -2.74576).abs() <= 1.0e-8);
    /// assert!((res6.re - -244.141).abs() <= 1.0e-3);
    /// assert!((res6.im - -378.594).abs() <= 1.0e-3);
    /// ```
    pub fn euler(n: usize) -> Self {

        // Initializing the vectors
        let mut coef: HashMap<i32, f64> = (0..=n).map(|x| (x as i32, 0.0)).collect();

        for k in 0..=n {
            let binom: f64 = basic::binomial(n, k) as f64;
            let f: f64 = Self::euler_number(k) / 2.0_f64.powi(k as i32);

            // Second loop from the (x - 1/2)^(n-k)
            for p in 0..=(n-k) {
                let pre: f64 = (-0.5_f64).powi(p as i32);
                let triangle_val: f64 = basic::binomial(n-k, p) as f64;
                let id: i32 = (n-k-p) as i32;

                if let Some(fact) = coef.get_mut(&id) {
                    *fact += binom * f * pre * triangle_val;
                }
            }
        }

        // Returning associated struct
        Self {
            coef,
            l: None,
            pre_f: None,
            compute_fn: Self::compute_base,
            compute_fnc: Self::compute_base_complex
        }
    }

    /// # Rising factorial polynomial
    /// 
    /// ## Definition
    /// The [rising factorial](https://en.wikipedia.org/wiki/Falling_and_rising_factorials) is
    /// defined as:
    /// $$
    /// x^{\overline{n}} = \prod_{k=0}^{n-1}(x+k)
    /// $$
    /// 
    /// The associated $n^{th}$ order polynomials can be computed using the [Stirling numbers](https://en.wikipedia.org/wiki/Stirling_numbers_of_the_first_kind),
    /// available in this crate as well.
    /// 
    /// ## Inputs
    /// - `n` the order of the polynomial
    /// 
    /// ## Example
    /// ```
    /// # use scilib::math::polynomial::Poly;
    /// let p5 = Poly::factorial_rising(5);
    /// let res = p5.compute(3.2);
    /// assert!((res - 3119.80032).abs() < 1e-8);
    /// ```
    pub fn factorial_rising(n: usize) -> Self {

        let mut coef: HashMap<i32, f64> = (0..=n).rev().map(|k| {
            (k as i32, Poly::stirling_number(n, k) as f64)
        }).collect();

        if n > 0 {
            coef.remove(&0);
        }

        Self {
            coef,
            l: None,
            pre_f: None,
            compute_fn: Self::compute_base,
            compute_fnc: Self::compute_base_complex
        }
    }
    
    /// # Falling factorial polynomial
    /// 
    /// ## Definition
    /// The [falling factorial](https://en.wikipedia.org/wiki/Falling_and_rising_factorials) is
    /// defined as:
    /// $$
    /// x^{\overline{n}} = \prod_{k=0}^{n-1}(x-k)
    /// $$
    /// 
    /// The associated $n^{th}$ order polynomials can be computed using the
    /// [signed Stirling numbers](https://en.wikipedia.org/wiki/Stirling_numbers_of_the_first_kind),
    /// available in this crate as well.
    /// 
    /// ## Inputs
    /// - `n` the order of the polynomial
    /// 
    /// ## Example
    /// ```
    /// # use scilib::math::polynomial::Poly;
    /// let p5 = Poly::factorial_falling(5);
    /// let res = p5.compute(3.2);
    /// assert!((res - -1.35168).abs() < 1e-8);
    /// ```
    pub fn factorial_falling(n: usize) -> Self {

        let mut coef: HashMap<i32, f64> = (0..=n).rev().map(|k| {
            (k as i32, Poly::stirling_number_signed(n, k) as f64)
        }).collect();

        if n > 0 {
            coef.remove(&0);
        }

        Self {
            coef,
            l: None,
            pre_f: None,
            compute_fn: Self::compute_base,
            compute_fnc: Self::compute_base_complex
        }
    }

    //////////////////////////////////////////////////
    // Polynomial operations

    /// # Polynomial derivation
    /// Computes the m derivative of the polynomial
    ///
    /// This is used to produce the $P_n^m(x)$ version of the polynomial.
    fn derive(&mut self, m: i32) {
        
        // Looping the derivation, m times
        for _ in 1..=m {

            let mut temp_c: HashMap<i32, f64> = HashMap::new();

            for (p, f) in &self.coef {
                match p {
                    0 => { },
                    _ => { temp_c.insert(p - 1, f * *p as f64); }
                }
            }
            self.coef = temp_c;
        }
    }

    /// # Computing for a real value
    /// 
    /// ## Definition
    /// Since there can be behaviors associated to some polynomials,
    /// we may need to call a specific computation function. To avoid hazards, we provide
    /// a single callable function that points to the correct computation function
    /// for a given polynomial.
    /// 
    /// ## Inputs
    /// - `x`: the value to evaluate (`x`: real)
    /// 
    /// Returns the value of the polynomial for a given $x$.
    pub fn compute(&self, x: f64) -> f64 {
        (self.compute_fn)(&self, x)
    }

    /// # Compute the polynomial for a real number
    fn compute_base(&self, x: f64) -> f64 {

        // Iterates through the values of the factors and powers
        self.coef.iter().fold(0.0, |res, (p, f)| res + f * x.powi(*p))
    }

    /// # Computation for Legendre
    fn compute_legendre(&self, x: f64) -> f64 {

        // Iterates through the values of the factors and powers
        let pre: f64 = self.pre_f.unwrap() * (1.0 - x.powi(2)).powf(self.l.unwrap() / 2.0);
        self.coef.iter().fold(0.0, |res, (p, f)| res + f * x.powi(*p)) * pre
    }
    /// # Computing for a complex value
    /// 
    /// ## Definition
    /// Since there can be behaviors associated to some polynomials,
    /// we may need to call a specific computation function. To avoid hazards, we provide
    /// a single callable function that points to the correct computation function
    /// for a given polynomial.
    /// 
    /// ## Inputs
    /// - `z`: the value to evaluate (`z`: complex)
    /// 
    /// Returns the value of the polynomial for a given $z$.
    pub fn compute_complex(&self, z: Complex64) -> Complex64 {
        (self.compute_fnc)(&self, z)
    }

    /// # Compute the polynomial for a complex number
    fn compute_base_complex(&self, z: Complex64) -> Complex64 {

        // Iterates through the values of the factors and powers
        self.coef.iter().fold(Complex64::default(), |res, (p, f)| res + f * z.powi(*p))
    }

    /// # Computation for Legendre
    fn compute_legendre_complex(&self, z: Complex64) -> Complex64 {

        // Iterates through the values of the factors and powers
        let pre: Complex64 = self.pre_f.unwrap() * (1.0 - z.powi(2)).powf(self.l.unwrap() / 2.0);
        self.coef.iter().fold(Complex64::default(), |res, (p, f)| res + f * z.powi(*p)) * pre
    }

    //////////////////////////////////////////////////
    // Specific methods to generate coefficients

    /// # Bernoulli number generator
    /// ## Definition
    /// Following the $B_m^{-}$ convention, we use the formula:
    /// $$
    /// B_m^{-} = \sum_{k=0}^{m}\sum_{v=0}^{k}(-1)^v \binom{k}{v}\frac{v^m}{k+1}
    /// $$
    /// 
    /// ## Inputs
    /// - `m`: the index of the number to compute
    /// 
    /// Returns the Associated Bernoulli number.
    /// 
    /// ## Example
    /// ```
    /// # use scilib::math::polynomial::Poly;
    /// let num_0: f64 = Poly::bernoulli_number(0);
    /// let num_1: f64 = Poly::bernoulli_number(1);
    /// let num_2: f64 = Poly::bernoulli_number(2);
    /// 
    /// assert!((num_0 - 1.0).abs() <= 1.0e-8);
    /// assert!((num_1 - -0.5).abs() <= 1.0e-8);
    /// assert!((num_2 - 1.0 / 6.0).abs() <= 1.0e-8);
    /// ```
    pub fn bernoulli_number(m: usize) -> f64 {

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

    /// # Euler number generator
    /// ## Definition
    /// We use the formula:
    /// $$
    /// E_k = \sum_{i=1}^{k}(-1)^{i}\frac{1}{2^i}\sum_{l=0}^{2i}(-1)^l\binom{2i}{l}(i-l)^{k}
    /// $$
    /// For every even $k$; else the result is 0.
    /// 
    /// ## Inputs
    /// - `m`: the index of the number to generate
    /// 
    /// Returns the associated Euler number.
    /// 
    /// ## Example
    /// ```
    /// # use scilib::math::polynomial::Poly;
    /// // We know that all odd Euler numbers are 0
    /// for odd in vec![1, 3, 5, 7, 9, 11] {
    ///     let ek: f64 = Poly::euler_number(odd);
    ///     assert_eq!(ek, 0.0);
    /// }
    /// 
    /// // The others oscillate between signs
    /// let e0: f64 = Poly::euler_number(0);
    /// let e2: f64 = Poly::euler_number(2);
    /// let e10: f64 = Poly::euler_number(10);
    /// 
    /// assert_eq!(e0, 1.0);
    /// assert_eq!(e2, -1.0);
    /// assert_eq!(e10, -50_521.0);
    /// ```
    pub fn euler_number(m: usize) -> f64 {

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

    /// # Stirling number generator
    /// 
    /// ## Definition
    /// The [Stirling numbers](https://en.wikipedia.org/wiki/Stirling_numbers_of_the_first_kind) are
    /// defined by the following recurrence:
    /// $$
    /// S(n, k) = (n - 1)S(n - 1, k) + S(n, k - 1)
    /// $$
    /// With the special conditions that:
    /// $$
    /// S(n, 0) = S(0, k) = 0,~ S(0, 0) = 1
    /// $$
    /// 
    /// ## Inputs
    /// - `n`: the order of the numbers
    /// - `k`: the position in the row
    /// 
    /// Returns the $S(n, k)$ Stirling number.
    /// 
    /// ## Example
    /// We can generate the fifth row of the Stirling numbers.
    /// 
    /// ```
    /// # use scilib::math::polynomial::Poly;
    /// let res: Vec<usize> = (0..=5).map(|k| Poly::stirling_number(5, k)).collect();
    /// let expected: Vec<usize> = vec![0, 24, 50, 35, 10, 1];
    /// assert_eq!(res, expected);
    /// ```
    pub fn stirling_number(n: usize, k: usize) -> usize {
        
        // We start with the special cases
        if n == 0 && k == 0 {
            1
        } else if n == 0 || k == 0 {
            0
        } else {
            (n - 1) * Poly::stirling_number(n - 1, k) + Poly::stirling_number(n - 1, k - 1)
        }
    }

    /// # Signed Stirling number generator
    /// 
    /// ## Definition
    /// The [signed Stirling numbers](https://en.wikipedia.org/wiki/Stirling_numbers_of_the_first_kind)
    /// have the same definition as the Stirling numbers (see definition of the function), but have
    /// a sign following the parity of the inputs:
    /// $$
    /// s(n, k) = (-1)^{n-k}S(n, k)
    /// $$
    /// 
    /// ## Inputs
    /// - `n`: the order of the numbers
    /// - `k`: the position in the row
    /// 
    /// Returns the $s(n, k)$ signed Stirling number.
    /// 
    /// ## Example
    /// We can generate the fourth row of the signed Stirling numbers.
    /// The values are similar to that of the Stirling numbers, but
    /// have some negatives.
    /// 
    /// ```
    /// # use scilib::math::polynomial::Poly;
    /// let res: Vec<i32> = (0..=4).map(|k| Poly::stirling_number_signed(4, k)).collect();
    /// let expected: Vec<i32> = vec![0, -6, 11, -6, 1];
    /// assert_eq!(res, expected);
    /// ```
    pub fn stirling_number_signed(n: usize, k: usize) -> i32 {
        // We start with the special cases
        let res = Self::stirling_number(n, k) as i32;

        (-1.0_f64).powi((n - k) as i32) as i32 * res
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Implementing operations

impl<T: Into<f64>> std::ops::Add<T> for Poly {
    type Output = Self;
    fn add(self, rhs: T) -> Self::Output {

        let mut coef = self.coef.clone();
        match coef.get_mut(&0) {
            Some(v) => *v += rhs.into(),
            None => { coef.insert(0, rhs.into()); }
        }

        Self {
            coef,
            ..self
        }
    }
}

impl<T: Into<f64>> std::ops::AddAssign<T> for Poly {
    fn add_assign(&mut self, rhs: T) {
        match self.coef.get_mut(&0) {
            Some(v) => *v += rhs.into(),
            None => { self.coef.insert(0, rhs.into()); }
        }
    }
}

impl<T: Into<f64>> std::ops::Sub<T> for Poly {
    type Output = Self;
    fn sub(self, rhs: T) -> Self::Output {

        let mut coef = self.coef.clone();
        match coef.get_mut(&0) {
            Some(v) => *v -= rhs.into(),
            None => { coef.insert(0, -rhs.into()); }
        }

        Self {
            coef,
            ..self
        }
    }
}

impl<T: Into<f64>> std::ops::SubAssign<T> for Poly {
    fn sub_assign(&mut self, rhs: T) {
        match self.coef.get_mut(&0) {
            Some(v) => *v -= rhs.into(),
            None => { self.coef.insert(0, -rhs.into()); }
        }
    }
}

impl<T: Into<f64>> std::ops::Mul<T> for Poly {
    type Output = Self;
    fn mul(self, rhs: T) -> Self::Output {

        let rhs_conv: f64 = rhs.into();
        let mut coef = self.coef.clone();

        for f in coef.values_mut() {
            *f *= rhs_conv;
        }
        
        Self {
            coef,
            ..self
        }
    }
}

impl<T: Into<f64>> std::ops::MulAssign<T> for Poly {
    fn mul_assign(&mut self, rhs: T) {
        let rhs_conv: f64 = rhs.into();
        for f in self.coef.values_mut() {
            *f *= rhs_conv;
        }
    }
}

impl<T: Into<f64>> std::ops::Div<T> for Poly {
    type Output = Self;
    fn div(self, rhs: T) -> Self::Output {

        let rhs_conv: f64 = rhs.into();
        let mut coef = self.coef.clone();

        for f in coef.values_mut() {
            *f /= rhs_conv;
        }
        
        Self {
            coef,
            ..self
        }
    }
}

impl<T: Into<f64>> std::ops::DivAssign<T> for Poly {
    fn div_assign(&mut self, rhs: T) {
        let rhs_conv: f64 = rhs.into();
        for f in self.coef.values_mut() {
            *f /= rhs_conv;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
