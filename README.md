# Scilib

> A rust crate for science

# Overview

This crate is designed to help any mathematical or scientific process for the Rust community. It compiles many useful concepts and items that are key in scientific applications, such as Bessel functions, statistical analysis, physical constants, etc...

The aim is to provide classical functions in pure Rust, for ease of operability.

---

## Work in progress

As of the creation of this readme, I am working on this project alone which means a few things:

1. The progression will be linked to my schedule
2. I will work firsts on concept with which I am familiar with
3. I am a self-taught developer, some solutions could be sub-optimal and thus improved

## What's coming?

The schedule of the development of the crate is not clear, as I am for now writing this as a side project. I plan on adding many useful functions from a physics point of view, but will expand as I go. For now, my objectives are:

- Astrophysics
- Thermodynamics
- Quantum mechanics (and statistical physics)
- Electromagnetism

And hopefully more when this is done (statistics for example)!

---

# Contents

## Complex numbers

This crate provides basic functionalities for complex numbers, mainly to support its other goals. The implementation uses `f64` for both the real and imaginary parts, to ensure precision in the computations.

Basic operations have been implemented to facilitate their use, and should be pretty easy to manipulate.

```rust
let c1 = Complex::from(2, 3.5);
let c2 = Complex::from(-1.2, 4) * 2;
println!("{}", c1 + c2);
```

***More functionalities are on their way, they will be added as they are needed for other domains.***

---

## Bessel functions

Essential in many maths and physics domain, bessel function are solutions of Bessel's differential equation ([Wiki page](https://en.wikipedia.org/wiki/Bessel_function)). This crate provides functions for both real and complex numbers, and for integer or real function order.

All functions are implemented:
- **J**: First kind
- **Y**: Second Kind
- **I**: Modified first kind
- **K**: Modified second kind
- **H1**: Hankel first kind
- **H2**: Hankel second kind

```rust
// All functions support complex numbers, and real orders
let res = bessel::jf(-1.2, 2.3);        // Computes -1.2 with order 2.3 in J
let res = bessel::y(3.5, 1);            // Y computes the limit for integer order
let res = bessel::hankel_first(2, -2)   // Hankel first kind
```

Values are compared to known results (thanks, [WolframAlpha](https://www.wolframalpha.com/)), and the results are within small margins of error.

***Documentation and tests are missing for I, K, H1 and H2, it's coming.***

---

## Typical polynomials

Useful polynomials will be implemented to facilitate their uses to everyone; as it stands, both the Legendre and Laguerre polynomials have been implemented.

```rust
// Legendre supports derivative (and negative m)
let leg = polynomial::Legendre::new(2, 1);  // l=2, m=1

// So does Laguerre
let lag = polynomial::Laguerre::new(3, -2); // l=3, m=-2
```

---

## Spherical harmonics

A first draft of a working solution for the spherical harmonics equation is implemented, and will be improved.

---
