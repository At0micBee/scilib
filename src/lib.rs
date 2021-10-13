//!
//! ## Introduction
//! 
//! > A rust crate for mathematics and science
//! 
//! This library aims at improving the state of the scientific viability of Rust. A number of key function
//! are missing that are hindering the use of the language in certain domains, despite the numerous
//! advantages of the language.
//! 
//! This crate will provide both simple, and advanced concept (from constants to spherical harmonics functions)
//! that will hopefully help many science developers.
//! 
//! > **Before you use**: this crate is currently a work in progress and is thus missing many features. I will do my best to ensure
//! both fast and correct computation, but it is evident that improvements could be found in the future. I will
//! first implement concept I am familiar with, and work towards other domains later on.
//! 
//! ## Contents
//! 
//! ### Useful mathematics function
//! 
//! The Rust library doesn't provide some functions that are quite common in scientific processes, and this crate attempts to provide as many as it can. Euler's Gamma and Beta function, Newton's binomial, factorial, the error functions (erf, erfc, erfi), ...
//! 
//! Some functions are still missing as of the writing of this document, but will be added later on.
//! 
//! ```rust
//! // These functions can be found in the math crate
//! use scilib::math::basic::*;
//! 
//! let g = gamma(3.2);
//! let b = beta(-1.2, 2.5);
//! 
//! // The erf function can compute Complex numbers (erfc, erfi as well)
//! let c = Complex::from(-0.1, 0.7);
//! let e = erf(c);
//! ```
//! ---
//! 
//! ### Coordinate systems
//! 
//! This crate provides functionalities for coordinate systems, such as Cartesian and Spherical, with many standard operations, and conversions.
//! 
//! ```rust
//! // They are found in the coordinate crate
//! use scilib::coordinate::*;
//! 
//! let c = cartesian::Cartesian::from(2.0, 1, 0.25);
//! let s = spherical::Spherical::from_degree(1.2, 30, 60.2);
//! ```
//! 
//! ---
//! 
//! ### Complex numbers
//! 
//! This crate provides basic functionalities for complex numbers, mainly to support its other goals. The implementation uses `f64` for both the real and imaginary parts, to ensure precision in the computations.
//! 
//! Basic operations have been implemented to facilitate their use, and should be pretty easy to manipulate.
//! 
//! ```rust
//! // They are found in the complex crate
//! use scilib::math::complex::Complex;
//! 
//! let c1 = Complex::from(2, 3.5);
//! let c2 = Complex::from(-1.2, 4) * 2;
//! println!("{}", c1 + c2);
//! ```
//! 
//! ***More functionalities are on their way, they will be added as they are needed for other domains.***
//! 
//! ---
//! 
//! ### Bessel functions
//! 
//! Essential in many maths and physics domain, bessel function are solutions of Bessel's differential equation ([Wiki page](https://en.wikipedia.org/wiki/Bessel_function)). This crate provides functions for both real and complex numbers, and for integer or real function order.
//! 
//! All functions are implemented:
//! - **J**: First kind
//! - **Y**: Second Kind
//! - **I**: Modified first kind
//! - **K**: Modified second kind
//! - **H1**: Hankel first kind
//! - **H2**: Hankel second kind
//! 
//! ```rust
//! // Found in the math crate
//! use scilib::math::bessel;
//! 
//! // All functions support complex numbers, and real orders
//! let res = bessel::jf(-1.2, 2.3);        // Computes -1.2 with order 2.3 in J
//! let res = bessel::y(3.5, 1);            // Y computes the limit for integer order
//! let res = bessel::hankel_first(2, -2)   // Hankel first kind
//! ```
//! 
//! Values are compared to known results (thanks, [WolframAlpha](https://www.wolframalpha.com/)), and the results are within small margins of error.
//! 
//! ---
//! 
//! ### Typical polynomials
//! 
//! Useful polynomials will be implemented to facilitate their uses to everyone; as it stands, both the Legendre (Plm(x)) and Laguerre (Llm(x)) polynomials have been implemented, where -l <= m <= l.
//! 
//! ```rust
//! // They are found in the polynomial crate
//! use scilib::math::polynomial;
//! 
//! // Legendre supports derivative (and negative m)
//! let leg = polynomial::Legendre::new(2, 1);  // l=2, m=1
//! 
//! // So does Laguerre
//! let lag = polynomial::Laguerre::new(3, -2); // l=3, m=-2
//! ```
//! 
//! ---
//! 
//! ### Quantum mechanics
//! 
//! The spherical harmonics Ylm(theta, phi) function has been added to the quantum section, and is valid for acoustics as well.
//! 
//! ```rust
//! // Found in the quantum crate
//! use scilib::quantum;
//! 
//! // Computing Ylm for l=3, m=1, theta = 0.2 and phi = -0.
//! 
//! let res = spherical_harmonics(3, 1, 0.2, -0.3);
//! ```
//! 
//! ---

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

pub mod astronomy;

pub mod coordinate;

pub mod constant;

pub mod math;

pub mod quantum;