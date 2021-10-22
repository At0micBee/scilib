//!
//! ## Introduction
//! 
//! > **A rust crate for mathematics and science**
//! 
//! This library aims at improving the state of the scientific viability of Rust. A number of key functions
//! are missing that are hindering the use of the language in certain domains, despite its other numerous advantages.
//! 
//! This crate will provide practical concepts (from constants to spherical harmonics functions)
//! that will hopefully help many science developers.
//! 
//! > **Before you use**: this crate is currently a work in progress and is thus missing many features. I will do my best to ensure
//! both fast and correct computation, but it is evident that improvements could be found in the future. I will
//! first implement concept I am familiar with, and work towards other domains later on. 
//! 
//! ## Contents
//! 
//! The scilib crate is sub-divided in themes, to simplify its use.
//! 
//! ### General purpose
//! 
//! - **Math**: Provides many base utilities, from complex numbers to bessel functions.
//! - **Coordinate**: Provides support for coordinate systems, and their respective operations.
//! - **Constant**: Contains many useful constants for physics
//! - **Signal**: Convolution and fast Fourier transform functions
//! 
//! ### Specific purpose
//! 
//! - **Astronomy**: Astronomical positions (***wip***)
//! - **Quantum**: Quantum mechanics functions (***wip***)
//! 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



pub mod astronomy;

pub mod constant;

pub mod coordinate;

pub mod math;

pub mod quantum;

pub mod range;

pub mod signal;
