//!
//! # Polynomials
//! 
//! Support to create and use Polynomials in Rust. They can be create either manually or by using
//! a function to generate one of many classical polynomial forms, such as Legendre, Laguerre, Hermite, ...
//! 
//! ## Definition
//! 
//! Polynomials are "simple" mathematical constructs that have the form:
//! $$
//! P_n(x) = \sum_{p=0}^{n} a_nx^n
//! $$
//! 
//! Where $x^n$ is the input value raised to the nth power, and the $a_n$ are associated coefficients to
//! each power. We can write the above example in an expanded form as:
//! $$
//! P_n(x) = a_0 + a_1x + a_2x^2 + ... + a_nx^n
//! $$
//! 
//! This crate allows for the creation and handling of polynomials of up to degree $n$, for example:
//! 
//! ```
//! use scilib::math::polynomial::Poly;
//! let p = Poly::from(&[(0, 1.0), (1, -1.0), (2, 2.0)]);
//! ```
//! 
//! creates the polynomial $1 - x + 2x^2$.
//! 
//! Basic operations are implemented, such as additions, subtraction, multiplication, both with numbers
//! and other polynomials:
//! 
//! ```
//! # use scilib::math::polynomial::Poly;
//! let p1 = Poly::from(&[(0, 1.0), (1, -1.0), (2, 2.0)]);
//! let p2 = Poly::from(&[(0, -2.0), (2, 1.2), (4, -0.2)]);
//! let p3 = p1 * 3.0;
//! let p4 = p2 / 2.0 + 1.5;
//! let mut res = p3 * p4;
//! res.derive(1);
//! ```
//! 
//! ## Implementation of named polynomials
//! 
//! To simplify the creation and use of typical polynomials, a variety of polynomials have been implemented.
//! So far, the list is:
//! 
//! - **Legendre**: `L(n,l)` generalized with with `n` positive integer and `l` positive or negative integer such that `-n < l < n`
//! - **Laguerre**: `L(n,l)` generalized with `n` positive integer and `l` a real number
//! - **Bernoulli**: `B(n)` with `n` positive integer
//! - **Euler**: `E(n)` with `n` positive integer
//! - **Bessel**: `y(n)` with `n` positive integer
//! - **Hermite**: `H(n)` with `n` positive integer
//! - **Rising factorial**: the polynomial associated to the rising factorial function, with `n` positive integer
//! - **Falling factorial**: the polynomial associated to the falling factorial function, with `n` positive integer
//! 
//! For example, to create the generalized Legendre Polynomial of degree 4, with associated factor 1:
//! 
//! ```
//! # use scilib::math::polynomial::Poly;
//! # use num_complex::Complex64;
//! let l = Poly::legendre(4, 1);   // And you're done!
//! 
//! // You can now use it to compute whatever you might want
//! let res = l.compute(2.7);
//! 
//! // Support for complex number using num_complex
//! let res_c = l.compute_complex(Complex64::from(1.0, 0.2));
//! ```
//! 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

use std::collections::HashMap;  // Rust Hashmap

use super::basic;               // Basic functions
use super::series;

use num_complex::Complex64;     // Using complex numbers from the num crate

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # Polynomial implementation
#[derive(Clone)]
pub struct Poly {
    coef: HashMap<i32, f64>,
    l: Option<f64>,
    compute_fn: fn(&Self, f64) -> f64,
    compute_fnc: fn(&Self, Complex64) -> Complex64
}

/// # Debug for polynomial
impl std::fmt::Debug for Poly {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        writeln!(f, "{:?}", self.coef)?;
        writeln!(f, "{:?}", self.l)?;

        Ok(())
    }
}

/// # Display for polynomials
impl std::fmt::Display for Poly {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        
        // We write the power and associated factors
        writeln!(f, "Power :: Factor")?;
        for (power, factor) in &self.coef {
            writeln!(f, "{:5} :: {}", power, factor)?;
        }

        Ok(())
    }
}

// # Default value for Poly
impl Default for Poly {
    fn default() -> Self {
        Self {
            coef: HashMap::new(),
            l: None,
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
    pub fn from(pow_fac: &[(usize, f64)]) -> Self {

        let coef: HashMap<i32, f64> = pow_fac.iter().map(|(p, f)| (*p as i32, *f)).collect();

        Self {
            coef,
            ..Self::default()
        }
    }

    //////////////////////////////////////////////////
    // Creating special polynomials

    /// # Legendre polynomials
    /// 
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
    /// To get the "regular" Legendre polynomial, simply set `l=0` when calling the function.
    /// 
    /// ## Inputs
    /// Produces the factors and powers for nth order polynomial, where
    /// - `n`: the order of the Legendre polynomial
    /// - `l`: the derivative order
    /// 
    /// /!\ IMPORTANT NOTE: when using a non-zero `l`, Legendre requires a secondary
    /// term to be multiplied to the polynomial during computation. This added term
    /// cannot be simplified and integrated to the base polynomial. If you need to 
    /// use the generalized version, be careful that arithmetic operations will
    /// be broken for it!
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
            compute_fn: Self::compute_legendre,
            compute_fnc: Self::compute_legendre_complex
        };

        // We derive the polynomial m times
        poly.derive(l.abs() as usize);

        // Returning the final polynomial
        poly * pre_f
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
            ..Self::default()
        }
    }

    /// # Bernoulli polynomials
    /// 
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
            ..Self::default()
        }
    }

    /// # Euler polynomials
    /// 
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
            ..Self::default()
        }
    }

    /// # Rising factorial polynomials
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
            ..Self::default()
        }
    }
    
    /// # Falling factorial polynomials
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
            ..Self::default()
        }
    }

    /// # Bessel Polynomials
    /// 
    /// ## Definition
    /// The [Bessel polynomials](https://en.wikipedia.org/wiki/Bessel_polynomials) are an orthogonal sequence of polynomials.
    /// They are defined as:
    /// $$
    /// y_n(x) = \sum_{k=0}^{n} x^k \frac{(n+k)!}{2^k k! (n-k)!}
    /// $$
    /// 
    /// The polynomials are related to the bessel functions through:
    /// $$
    /// y_n(x) = \sqrt{\frac{2}{x\pi}}e^{1/x}K_{n+\frac{1}{2}}\left(\frac{1}{x}\right)
    /// $$
    /// where $K$ is a modified Bessel function of the second kind, implemented under the `math::bessel` module.
    /// 
    /// ## Inputs
    /// - `n` the order of the polynomial
    /// 
    /// ## Example
    /// ```
    /// # use scilib::math::polynomial::Poly;
    /// # use scilib::math::bessel;
    /// use std::f64::consts::FRAC_PI_2;
    /// let x: f64 = 0.28;
    /// let order: usize = 5;
    /// 
    /// let p = Poly::bessel(order);    // n = 5
    /// let res = p.compute(x);
    /// assert!((res - 30.086718976000007).abs() < 1e-8);
    /// 
    /// // We can also check that the results are coherent with the values in math::bessel
    /// let fact: f64 = (1.0 / (x * FRAC_PI_2)).sqrt() * (1.0 / x).exp();
    /// let res_b = fact * bessel::k(1.0 / x, order as f64 + 0.5);
    /// assert!((res - res_b.re).abs() <= 1e-3);
    /// ```
    pub fn bessel(n: usize) -> Self {

        let mut coef: HashMap<i32, f64> = HashMap::new();

        // Setting up variables once
        let mut kf: usize = 1;                      // k!
        let mut twos: usize = 1;                    // 2^k
        let mut top: usize = basic::factorial(n);   // (n + k)!, computing the factorial only once
        let mut bot: usize = top;                   // (n - k)!
        let mut c: f64;                             // The result for the coef
        coef.insert(0, 1.0);                        // Value for 0 is always 1

        // Starting iterations at 1
        for k in 1..=n {

            kf *= k;
            top *= n + k;
            twos *= 2;
            bot /= n - k + 1;
            c = top as f64 / (twos * bot * kf) as f64;

            coef.insert(k as i32, c);
        }

        Self {
            coef,
            ..Self::default()
        }
    }

    /// # Hermite polynomials
    /// 
    /// ## Definitions
    /// The [Hermite polynomials](https://en.wikipedia.org/wiki/Hermite_polynomials) are a sequence
    /// of orthogonal polynomials, widely used in a variety of applications.
    /// This functions generates the so called "physicist" version of the polynomials, defined as:
    /// $$
    /// H_n(x) = (-1)^n e^{x^2} \frac{de^{-x^2}}{dx^n}
    /// $$
    /// 
    /// The polynomials coefficients are computed using the following generating function:
    /// $$
    /// H_n(x) = n! \sum_{m=0}^{\lfloor\frac{n}{2}\rfloor}  \frac{(-1)^m 2^{n-2m}}{m!(n-2m)!} x^{n-2m}
    /// $$
    /// 
    /// ## Inputs
    /// 
    /// - `n` the order of the polynomial
    /// 
    /// ## Example
    /// ```
    /// # use scilib::math::polynomial::Poly;
    /// let h = Poly::hermite(7);
    /// let res = h.compute(1.2);
    /// assert!((res - 904.4250624000001).abs() < 1e-8);
    /// ```
    pub fn hermite(n: usize) -> Self {

        let mut coef: HashMap<i32, f64> = HashMap::new();

        let nf: usize = basic::factorial(n);
        let mut n2m: usize;
        let mut sign: isize = 1;
        let mut twos: usize = 2_usize.pow(n as u32);
        let mut mf: usize = 1;
        let mut c: f64;
        coef.insert(n as i32, twos as f64);
        
        for m in 1..=(n / 2) {

            sign *= -1;
            twos = 2_usize.pow((n - 2 * m) as u32);
            mf *= m;
            n2m = basic::factorial(n - 2 * m);
            c = sign as f64 * (nf * twos) as f64 / (mf * n2m) as f64;

            coef.insert((n - 2 * m) as i32, c);
        }

        Self {
            coef,
            ..Self::default()
        }
    }

    //////////////////////////////////////////////////
    // Polynomial operations

    /// # Polynomial derivation
    /// 
    /// ## Definition
    /// Computes the mth derivative of the polynomial.
    /// 
    /// ## Inputs
    /// - `m` the derivative order to compute
    ///
    /// ## Example
    /// ```
    /// # use scilib::math::polynomial::Poly;
    /// let mut p = Poly::from(&[(0, 1.0), (1, -1.0), (2, 2.0)]);
    /// let expected = Poly::from(&[(0, -1.0), (1, 4.0)]);
    /// p.derive(1);
    /// assert_eq!(p, expected);
    /// ```
    pub fn derive(&mut self, m: usize) {
        
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

    /// # Polynomial integration
    /// 
    /// ## Definition
    /// Computes the mth integration of the polynomial.
    /// 
    /// ## Inputs
    /// - `m` the integration order to compute
    ///
    /// ## Example
    /// ```
    /// # use scilib::math::polynomial::Poly;
    /// let mut p = Poly::from(&[(0, 1.0), (1, -1.0), (2, 2.0)]);
    /// let expected = Poly::from(&[(0, 2.2), (1, 1.0), (2, -0.5), (3, 2.0 / 3.0)]);
    /// p.integrate(1, &[2.2]);
    /// assert_eq!(p, expected);
    /// ```
    pub fn integrate(&mut self, m: usize, coef: &[f64]) {

        // We check if coefficients are passed
        let n_coef: Vec<f64> = if coef.len() > 0 {
            assert_eq!(coef.len(), m);  // If yes we check if there is the right number
            coef.into()                 // And we return the values
        } else {
            vec![1.0; m]                // If not passed then ones are used
        };
        
        // We loop to integrate
        for n_f in n_coef {

            let mut temp_c: HashMap<i32, f64> = HashMap::new();

            for (p, f) in &self.coef {
                temp_c.insert(p + 1, f / (p + 1) as f64);
            }

            temp_c.insert(0, n_f);
            
            self.coef = temp_c;
        }
    }

    /// # Raising polynomial to an integer power
    /// 
    /// ## Definition
    /// We compute the new polynomial that corresponds to the power it is raised to. 
    /// 
    /// ## Inputs
    /// - `n`: the power to raise the polynomial to
    /// 
    /// Returns: a new instance of `Poly`, which is the nth power of the current one.
    /// 
    /// ## Example
    /// ```
    /// # use scilib::math::polynomial::Poly;
    /// let p = Poly::from(&[(0, 1.0), (2, 2.2), (3, -1.2)]);
    /// let n = p.pow(3);
    /// let n_coefs = n.get_coefs();
    /// 
    /// let theory = vec![
    ///     (0, 1.0), (2, 6.6), (3, -3.6), (4, 14.52), (5, -15.84),
    ///     (6, 14.968), (7, -17.424), (8, 9.504), (9, -1.728)
    /// ];
    /// let expected = Poly::from(&theory);
    /// let expected_coefs = expected.get_coefs();
    /// 
    /// for (p, f) in &n_coefs {
    ///     match expected_coefs.get(p) {
    ///         Some(ef) => { assert!((f - ef).abs() < 1e-12) },
    ///         None => { assert!(false) }
    ///     }
    /// }
    /// 
    /// ```
    pub fn pow(&self, n: usize) -> Poly {

        if n == 0 {
            return Poly::from(&[(0, 1.0)]); // If n = 0
        }
        
        let mut res = self.clone();         // For n = 1
        
        for _ in 2..=n {                    // For n >= 2
            res *= self.clone();
        }

        res
    }

    /// # Getting the coefficients of a polynomial
    /// 
    /// Simply returns the HashMap of the power and associated coefficients.
    /// 
    /// ```
    /// # use std::collections::HashMap;
    /// # use scilib::math::polynomial::Poly;
    /// let setup = vec![(0, 1.0), (2, 2.2), (3, -1.2)];
    /// let p = Poly::from(&setup);
    /// let expected: HashMap<i32, f64> = setup.iter().map(|(p, f)| (*p as i32, *f)).collect();
    /// assert_eq!(p.get_coefs(), expected);
    /// ```
    pub fn get_coefs(&self) -> HashMap<i32, f64> {
        self.coef.clone()
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
        let pre: f64 = (1.0 - x.powi(2)).powf(self.l.unwrap() / 2.0);
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
        let pre: Complex64 = (1.0 - z.powi(2)).powf(self.l.unwrap() / 2.0);
        self.coef.iter().fold(Complex64::default(), |res, (p, f)| res + f * z.powi(*p)) * pre
    }

    /// # Polynomial orders
    /// 
    /// ## Definition
    /// Extracts the maximum order of the polynomial. Used built-in series function.
    /// 
    /// ## Example
    /// ```
    /// # use scilib::math::polynomial::Poly;
    /// let p72 = Poly::factorial_rising(7);    // n=7
    /// assert_eq!(p72.get_order(), 7);
    /// ```
    pub fn get_order(&self) -> i32 {
        // Collecting the powers
        let p: Vec<i32> = self.coef.keys().map(|x| *x).collect();

        // Computing the max in the slice
        series::max_slice(&p)
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

impl std::cmp::PartialEq for Poly {
    fn eq(&self, other: &Self) -> bool {

        let c: bool = self.coef == other.coef;
        let l: bool = self.l == other.l;

        c & l
    }
}

/// # Addition
/// 
/// ```
/// # use scilib::math::polynomial::Poly;
/// let mut p = Poly::from(&[(1, 2.0)]);
/// let mut e = Poly::from(&[(1, 2.0), (0, 5.0)]);
/// let k = p + 5;
/// assert_eq!(k, e);
/// ```
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

/// # Addition of Poly
/// 
/// ```
/// # use scilib::math::polynomial::Poly;
/// let p1 = Poly::from(&[(1, 2.5), (2, -1.0), (3, 1.2)]);
/// let p2 = Poly::from(&[(0, 10.0), (1, -5.0), (3, 1.2)]);
/// let res = p1 + p2;
/// let expected = Poly::from(&[(0, 10.0), (1, -2.5), (2, -1.0), (3, 2.4)]);
/// assert_eq!(res, expected);
/// ```
impl std::ops::Add<Self> for Poly {
    type Output = Self;
    fn add(self, rhs: Poly) -> Self::Output {

        let mut coef = self.coef.clone();
        for (c, f) in &rhs.coef {
            match coef.get_mut(c) {
                Some(v) => *v += f,
                None => { coef.insert(*c, *f); }
            }
        }
        
        Poly {
            coef,
            ..self
        }
    }    
}

/// # Addition assigned
/// 
/// ```
/// # use scilib::math::polynomial::Poly;
/// let mut p = Poly::from(&[(1, 2.0)]);
/// let mut e = Poly::from(&[(1, 2.0), (0, 5.0)]);
/// p += 5;
/// assert_eq!(p, e);
/// ```
impl<T: Into<f64>> std::ops::AddAssign<T> for Poly {
    fn add_assign(&mut self, rhs: T) {
        match self.coef.get_mut(&0) {
            Some(v) => *v += rhs.into(),
            None => { self.coef.insert(0, rhs.into()); }
        }
    }
}

/// # Addition assigned of Poly
/// 
/// ```
/// # use scilib::math::polynomial::Poly;
/// let mut p1 = Poly::from(&[(1, 2.5), (2, -1.0), (3, 1.2)]);
/// let p2 = Poly::from(&[(0, 10.0), (1, -5.0), (3, 1.2)]);
/// p1 += p2;
/// let expected = Poly::from(&[(0, 10.0), (1, -2.5), (2, -1.0), (3, 2.4)]);
/// assert_eq!(p1, expected);
/// ```
impl std::ops::AddAssign<Self> for Poly {
    fn add_assign(&mut self, rhs: Self) {
        for (c, f) in &rhs.coef {
            match self.coef.get_mut(c) {
                Some(v) => *v += f,
                None => { self.coef.insert(*c, *f); }
            }
        }
    }
}

/// # Subtraction
/// 
/// ```
/// # use scilib::math::polynomial::Poly;
/// let mut p = Poly::from(&[(1, 2.0)]);
/// let mut e = Poly::from(&[(1, 2.0), (0, -5.0)]);
/// let k = p - 5;
/// assert_eq!(k, e);
/// ```
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

/// # Subtraction of Poly
/// 
/// ```
/// # use scilib::math::polynomial::Poly;
/// let p1 = Poly::from(&[(1, 2.5), (2, -1.0), (3, 1.2)]);
/// let p2 = Poly::from(&[(0, 10.0), (1, -5.0), (3, 1.2)]);
/// let res = p1 - p2;
/// let expected = Poly::from(&[(0, -10.0), (1, 7.5), (2, -1.0), (3, 0.0)]);
/// assert_eq!(res, expected);
/// ```
impl std::ops::Sub<Self> for Poly {
    type Output = Self;
    fn sub(self, rhs: Poly) -> Self::Output {

        let mut coef = self.coef.clone();
        for (c, f) in &rhs.coef {
            match coef.get_mut(c) {
                Some(v) => *v -= f,
                None => { coef.insert(*c, -*f); }
            }
        }
        
        Poly {
            coef,
            ..self
        }
    }    
}

/// # Subtraction assigned
/// 
/// ```
/// # use scilib::math::polynomial::Poly;
/// let mut p = Poly::from(&[(1, 2.0)]);
/// let mut e = Poly::from(&[(1, 2.0), (0, -5.0)]);
/// p -= 5;
/// assert_eq!(p, e);
/// ```
impl<T: Into<f64>> std::ops::SubAssign<T> for Poly {
    fn sub_assign(&mut self, rhs: T) {
        match self.coef.get_mut(&0) {
            Some(v) => *v -= rhs.into(),
            None => { self.coef.insert(0, -rhs.into()); }
        }
    }
}

/// # Subtraction assigned of Poly
/// 
/// ```
/// # use scilib::math::polynomial::Poly;
/// let mut p1 = Poly::from(&[(1, 2.5), (2, -1.0), (3, 1.2)]);
/// let p2 = Poly::from(&[(0, 10.0), (1, -5.0), (3, 1.2)]);
/// p1 -= p2;
/// let expected = Poly::from(&[(0, -10.0), (1, 7.5), (2, -1.0), (3, 0.0)]);
/// assert_eq!(p1, expected);
/// ```
impl std::ops::SubAssign<Self> for Poly {
    fn sub_assign(&mut self, rhs: Self) {
        for (c, f) in &rhs.coef {
            match self.coef.get_mut(c) {
                Some(v) => *v -= f,
                None => { self.coef.insert(*c, -*f); }
            }
        }
    }
}

/// # Multiplication
/// 
/// ```
/// # use scilib::math::polynomial::Poly;
/// let mut p = Poly::from(&[(1, 2.0), (2, -1.0)]);
/// let mut e = Poly::from(&[(1, 6.0), (2, -3.0)]);
/// let k = p * 3;
/// assert_eq!(k, e);
/// ```
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

/// # Multiplication of Poly
/// 
/// ```
/// # use scilib::math::polynomial::Poly;
/// let p1 = Poly::from(&[(0, 1.0), (1, 2.5), (3, -1.2)]);
/// let p2 = Poly::from(&[(0, 1.0), (2, -1.0)]);
/// let res = p1 * p2;
/// let expected = Poly::from(&[(0, 1.0), (1, 2.5), (2, -1.0), (3, -3.7), (5, 1.2)]);
/// assert_eq!(res, expected);
/// ```
impl std::ops::Mul<Self> for Poly {
    type Output = Self;
    fn mul(self, rhs: Poly) -> Self::Output {

        let mut n_p: i32;
        let mut n_f: f64;

        let mut coef: HashMap<i32, f64> = HashMap::new();

        for (p, f) in self.coef.iter() {

            for (rhs_p, rhs_f) in &rhs.coef {
                n_p = p + rhs_p;
                n_f = f * rhs_f;

                match coef.get_mut(&n_p) {
                    Some(v) => *v += n_f,
                    None => { coef.insert(n_p, n_f); }
                }
            }
        };
        
        Poly {
            coef,
            ..self
        }
    }    
}

/// # Multiplication assigned
/// 
/// ```
/// # use scilib::math::polynomial::Poly;
/// let mut p = Poly::from(&[(1, 2.0), (2, -1.0)]);
/// let mut e = Poly::from(&[(1, 6.0), (2, -3.0)]);
/// p *= 3;
/// assert_eq!(p, e);
/// ```
impl<T: Into<f64>> std::ops::MulAssign<T> for Poly {
    fn mul_assign(&mut self, rhs: T) {
        let rhs_conv: f64 = rhs.into();
        for f in self.coef.values_mut() {
            *f *= rhs_conv;
        }
    }
}

/// # Multiplication assigned of Poly
/// 
/// ```
/// # use scilib::math::polynomial::Poly;
/// let mut p1 = Poly::from(&[(0, 1.0), (1, 2.5), (3, -1.2)]);
/// let p2 = Poly::from(&[(0, 1.0), (2, -1.0)]);
/// p1 *= p2;
/// let expected = Poly::from(&[(0, 1.0), (1, 2.5), (2, -1.0), (3, -3.7), (5, 1.2)]);
/// assert_eq!(p1, expected);
/// ```
impl std::ops::MulAssign<Self> for Poly {
    fn mul_assign(&mut self, rhs: Self) {
        let mut n_p: i32;
        let mut n_f: f64;

        let mut coef: HashMap<i32, f64> = HashMap::new();

        for (p, f) in self.coef.iter() {

            for (rhs_p, rhs_f) in &rhs.coef {
                n_p = p + rhs_p;
                n_f = f * rhs_f;

                match coef.get_mut(&n_p) {
                    Some(v) => *v += n_f,
                    None => { coef.insert(n_p, n_f); }
                }
            }
        };

        self.coef = coef;
    }
}

/// # Division
/// 
/// ```
/// # use scilib::math::polynomial::Poly;
/// let mut p = Poly::from(&[(1, 2.0), (2, -1.0)]);
/// let mut e = Poly::from(&[(1, 1.0), (2, -0.5)]);
/// let k = p / 2;
/// assert_eq!(k, e);
/// ```
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

/// # Division assigned
/// 
/// ```
/// # use scilib::math::polynomial::Poly;
/// let mut p = Poly::from(&[(1, 2.0), (2, -1.0)]);
/// let mut e = Poly::from(&[(1, 1.0), (2, -0.5)]);
/// p /= 2;
/// assert_eq!(p, e);
/// ```
impl<T: Into<f64>> std::ops::DivAssign<T> for Poly {
    fn div_assign(&mut self, rhs: T) {
        let rhs_conv: f64 = rhs.into();
        for f in self.coef.values_mut() {
            *f /= rhs_conv;
        }
    }
}

/// # Negation
/// 
/// ```
/// # use scilib::math::polynomial::Poly;
/// let mut p = Poly::from(&[(1, 2.0), (2, -1.0)]);
/// let mut e = Poly::from(&[(1, -2.0), (2, 1.0)]);
/// let k = -p;
/// assert_eq!(k, e);
/// ```
impl std::ops::Neg for Poly {
    type Output = Self;
    fn neg(self) -> Self::Output {
        
        self * -1
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
