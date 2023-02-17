
![](https://raw.githubusercontent.com/At0micBee/scilib/master/branding/Scilib.png)

---

# Overview

![](https://img.shields.io/docsrs/scilib?label=Tests&style=flat-square)
![](https://img.shields.io/crates/v/scilib?style=flat-square)
![](https://img.shields.io/crates/l/scilib?style=flat-square)

This crate is designed to help any mathematical or scientific processes for the Rust community. It compiles many useful concepts and items that are key in many scientific domains. The aim of the crate is to provide these functions in pure Rust, to avoid any dependencies.

---

# Contents

## Constants

Many useful constants have been added, comprising many different fields, from astrophysics to quantum mechanics, but also mathematics, thermodynamics, electromagnetism, etc... They're listed in the `constant` module. Note that all constants are provided with a link to the source.

```rust
use scilib::constant;

println!("{}", constant::SUN_RADIUS);   // Solar radius
println!("{}", constant::H_BAR);        // H bar
println!("{}", constant::K_B);          // Boltzmann constant
println!("{}", constant::BOHR_MAG);     // Bohr magneton
// And many more...
```

---

## Useful mathematical functions

The Rust library doesn't provide some functions that are quite common in scientific processes, and this crate attempts to provide as many as it can. Euler's Gamma and Beta function, Newton's binomial, factorial, the error functions (erf, erfc, erfi), ...

```rust
// These functions can be found in the math crate
use scilib::math::basic::*;

let g = gamma(3.2);
let b = beta(-1.2, 2.5);

// The erf function can compute Complex numbers (erfc, erfi as well)
let c = Complex64::new(-0.1, 0.7);
let e = erf(c);
```

---

## Bessel functions

Essential in many maths and physics domain, the **[bessel functions](https://en.wikipedia.org/wiki/Bessel_function)** are solutions of Bessel's differential equation. This crate provides functions for both real and complex numbers, and for integer or real function order. It covers standard Bessel functions, the spherical Bessel functions, and the Riccati-Bessel functions.

All functions are implemented:
- **J**: First kind
- **Y**: Second Kind
- **I**: Modified first kind
- **K**: Modified second kind
- **H1**: Hankel first kind
- **H2**: Hankel second kind
- **j**: Spherical first kind
- **y**: Spherical second kind
- **h1**: Spherical hankel first kind
- **h2**: Spherical hankel second kind
- **S**: Riccati-Bessel first kind
- **C**: Riccati-Bessel Second kind
- **Xi**: Riccati-Bessel Third kind
- **Zeta**: Riccati-Bessel Fourth kind

```rust
// Found in the math crate
use scilib::math::bessel;

// All functions support complex numbers, and real orders
let res_j = bessel::jf(-1.2, -2.3);             // J function; works for any input and order
let res_y = bessel::y(3.5, 1);                  // Y function; computes the limit for integer order
let res_i = bessel::i(7.2, 2.25);               // I function; similar to J
let res_k = bessel::k(-1.1, 0.5);               // K function; computes the limit for integer order
let res_1 = bessel::hankel_first(2, -2);        // Hankel first kind
let res_2 = bessel::hankel_second(1, -1.32);    // Hankel first kind
let res_sj = bessel::sj(4.4, 2);                // Spherical first kind
let res_sy = bessel::sy(-1.54, 3);              // Spherical second kind
let res_s1 = bessel::sh_first(2.11, 4);         // Spherical hankel first kind
let res_s2 = bessel::sh_second(0.253, 0);       // Spherical hankel second kind
```

---

## Typical polynomials

A dedicated method for polynomial is implemented in the module `math::polynomial` as `Poly`.
Many useful polynomials have also been implemented.

- **Legendre**: `L(n,l)` generalized with with `n` positive integer and `l` positive or negative integer such that `-n < l < n`
- **Laguerre**: `L(n,l)` generalized with `n` positive integer and `l` a real number
- **Bernoulli**: `B(n)` with `n` positive integer
- **Euler**: `E(n)` with `n` positive integer
- **Bessel**: `y(n)` with `n` positive integer
- **Rising factorial**: the polynomial associated to the rising factorial function, with `n` positive integer
- **Falling factorial**: the polynomial associated to the falling factorial function, with `n` positive integer

```rust
// They are found in the polynomial crate
use scilib::math::polynomial::Poly;

let p = Poly::from([(2, 1.0), (1, 2.0), (0, -1.0)]) // x² + 2x - 1

let leg = Poly::legendre(2, -1);                    // n=2, l=-1
let lag = Poly::laguerre(3, 2.78);                  // n=3, l=2.78
let ber = Poly::bernoulli(3);                       // n=3
let eul = Poly::euler(5);                           // n=5
let rf = Poly::factorial_rising(4);                 // n=4
let ff = Poly::factorial_falling(3);                // n=3
```

---

## Coordinate systems

This crate provides functionalities for coordinate systems, such as Cartesian and Spherical, with many standard operations and conversions.

```rust
// They are found in the coordinate crate
use scilib::coordinate::*;

let car = cartesian::Cartesian::from(2.0, 1, 0.25);
let sph = spherical::Spherical::from_degree(1.2, 30, 60.2);
let cyl = spherical::Cylindrical::from_degree(1.2, 30, -2.55);
```

---

## Signal functions

Support to conduct both fast Fourier transform (`fft`) and the inverse fast Fourier transform (`ifft`) is available. Computations are
done using [Bluestein's algorithm](https://en.wikipedia.org/wiki/Chirp_Z-transform#Bluestein.27s_algorithm). Convolution is also possible,
with any two vector sizes.

```rust
// Found in the fourier crate
use scilib::signal::*

// Computing values of the sinus
let r = range::linear(0.0, 10.0, 15);
let s: Vec<Complex64> = r.iter().map(|val| val.sin()).collect();

let res = fft(&s);
let res2 = ifft(&res);
let res3 = convolve(&r, &s);
```

---

## Astronomy and astrophysics

We provide practical functions for astronomy and astrophysics applications, from a Radec coordinate system to equilibrium temperature computation and a magnitude calculator.

```rust
// Found in the astronomy crate
use scilib::astronomy::*;
use scilib::constant as cst;

// Creating a Radec system
let coord: Radec = Radec::from_degree(32, 21.22534);

// And other practical function
let mag = apparent_mag(cst::SUN_L, cst::LY);            // Apparent mag of the Sun at 1 ly
let hill = hill_radius(mass, mass_star, distance, e);   // Hill radius
let b = impact_parameter(a, r_star, i, e, w);           // Transit impact parameter
```

---

## Quantum mechanics

Both the radial wave function Rnl(r) and the spherical harmonics Ylm(theta, phi) have been added to the quantum section. The Ylm is also valid for acoustics as well.

```rust
// Found in the quantum crate
use scilib::quantum::*;

// Computing Ylm for l=3, m=1, theta = 0.2 and phi = -0.3
let sph = spherical_harmonics(3, 1, 0.2, -0.3);

// Computing the Rnl for n=4, l=2
let rad = radial_wavefunction(4, 2, 1.3e-12);
```

---
