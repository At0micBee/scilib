//!
//! # Fourier transform algorithms
//! 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

use std::f64::consts::PI;

use crate::math::complex::Complex;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # Convolution
/// 
/// Computes the convolution of two vectors, including the edges.
/// 
/// ```
/// # use scilib::signal::convolve;
/// // Creating two vectors to convolve
/// let a1: Vec<f64> = vec![2.8, 2.5, 1.0, 0.5, 3.2, 0.25];
/// let a2: Vec<f64> = vec![0.25, 1.0, 0.5];
/// let res = convolve(&a1, &a2);
/// let expected: Vec<f64> = vec![0.7, 3.425, 4.15, 2.375, 1.8, 3.5125, 1.85, 0.125];
/// 
/// assert_eq!(res, expected);
/// ```
pub fn convolve<T>(a_i: &[T], b_i: &[T]) -> Vec<T>
where T: std::ops::Mul<Output = T> + std::ops::AddAssign + Default + Copy {

    // We check which box is the smallest
    let (a, b): (&[T], &[T]) = match a_i.len() < b_i.len() {
        true => (b_i, a_i),
        false => (a_i, b_i)
    };

    // We initialize our variables
    let mut sum: T;
    let l_a: usize = a.len();
    let l_b: usize = b.len();
    let mut res: Vec<T> = Vec::with_capacity(l_a + 2 * l_b - 2);    // Complete convolution length

    // Box b starts to slide over a
    for n in 0..(l_b - 1) {
        sum = T::default();
        for m in 0..=n {
            sum += a[n-m] * b[m];

        }
        res.push(sum);
    }

    // Box b covers a completely
    for n in (l_b - 1)..=(l_a - 1) {
        sum = T::default();
        for m in 0..l_b {
            sum += a[n-m] * b[m];
        }
        res.push(sum);
    }

    // Box b exits the slide over a
    for n in l_a..(l_a + l_b - 1) {
        sum = T::default();
        for m in (n - l_a + 1)..l_b {
            sum += a[n-m] * b[m];
        }
        res.push(sum);
    }

    res
}

/// # Exact convolution
/// 
/// Computes the convolution of two vectors, only for points where both vectors
/// are completely covered.
/// 
/// ```
/// # use scilib::signal::convolve_full;
/// // Creating two vectors to convolve
/// let a1: Vec<f64> = vec![2.8, 2.5, 1.0, 12.2];
/// let a2: Vec<f64> = vec![0.25, 1.0, 23.2, 0.25, 1.3, 2.1];
/// let res = convolve_full(&a1, &a2);
/// let expected: Vec<f64> = vec![62.75, 39.665, 292.42];
/// 
/// for (e, c) in expected.iter().zip(&res) {
///     assert!((e - c).abs() < 1.0e-10);
/// }
/// ```
pub fn convolve_full<T>(a_i: &[T], b_i: &[T]) -> Vec<T>
where T: std::ops::Mul<Output = T> + std::ops::AddAssign + Default + Copy {

    // We check which box is the smallest
    let (a, b): (&[T], &[T]) = match a_i.len() < b_i.len() {
        true => (b_i, a_i),
        false => (a_i, b_i)
    };

    // We initialize our variables
    let mut sum: T;
    let l_a: usize = a.len();
    let l_b: usize = b.len();
    let mut res: Vec<T> = Vec::with_capacity(l_a - l_b + 1);    // Exact convolution length

    // Box b covers a completely only
    for n in (l_b - 1)..=(l_a - 1) {
        sum = T::default();
        for m in 0..l_b {
            sum += a[n-m] * b[m];
        }
        res.push(sum);
    }

    res
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # Fast Fourier transform
/// 
/// Computes the FFT for a one-dimensional array, based on the discrete fourier transform.
/// This function accepts complex input. The FFT is computed using
/// [Bluestein's algorithm](https://en.wikipedia.org/wiki/Chirp_Z-transform#Bluestein.27s_algorithm).
/// 
/// ```
/// # use scilib::range;
/// # use scilib::math::complex::Complex;
/// # use scilib::signal::fft;
/// // We create a Vec with the value sin(v) + cos(v)i
/// let r = range::linear(0.0, 10.0, 15);
/// let s: Vec<Complex> = r.iter().map(|val| (Complex::i() * val.cos()) + val.sin()).collect();
/// let res = fft(&s);
/// 
/// // Checking for some values...
/// assert!((res[0].re - 2.19227873).abs() < 1.0e-8 && (res[0].im - -0.64850436).abs() < 1.0e-8);
/// assert!((res[4].re - 0.73245756).abs() < 1.0e-8 && (res[4].im - 0.44922173).abs() < 1.0e-8);
/// assert!((res[9].re - -0.02709553).abs() < 1.0e-8 && (res[9].im - 1.02037473).abs() < 1.0e-8);
/// assert!((res[14].re - 4.77371673).abs() < 1.0e-8 && (res[14].im - -2.58964065).abs() < 1.0e-8);
/// ```
pub fn fft<T>(data: &[T]) -> Vec<Complex>
where T: Into<Complex> + Copy {

    let length: isize = data.len() as isize;

    // The pre-factor in the exponential terms
    let ipi_n: Complex = Complex::i() * PI / length as f64;

    let mut res: Vec<Complex> = Vec::with_capacity(data.len());

    let mut an: Complex;
    let mut bkn: Complex;
    let mut bk_star: Complex;
    let mut sum_k: Complex;

    for k in 0..length {
        sum_k = Complex::default();
        bk_star = (k.pow(2) as f64 * ipi_n).exp().conjugate();

        for (n, val) in data.iter().enumerate() {
            an = (- (n.pow(2) as f64) * ipi_n).exp() * (*val).into();
            bkn = ((k - n as isize).pow(2) as f64 * ipi_n).exp();
            sum_k += an * bkn;
        }

        res.push(bk_star * sum_k);
    }

    res
}

/// # Inverse fast Fourier transform
/// 
/// Computes the IFFT for a one-dimensional array, based on the discrete fourier transform.
/// This function accepts complex input. The FFT is computed using
/// [Bluestein's algorithm](https://en.wikipedia.org/wiki/Chirp_Z-transform#Bluestein.27s_algorithm).
/// 
/// This function yields `v = ifft(fft(v))`, within numerical errors.
/// 
/// ```
/// # use scilib::range;
/// # use scilib::math::complex::Complex;
/// # use scilib::signal::{ fft, ifft };
/// // We create a Vec with the value sin(v) + cos(v)i
/// let r = range::linear(0.0, 10.0, 15);
/// let s: Vec<f64> = r.iter().map(|val| val.sin()).collect();
/// let res = ifft(&fft(&s));
/// 
/// for (ori, comp) in s.iter().zip(&res) {
///     assert!((ori - comp.re).abs() < 1.0e-14 && comp.im < 1.0e-14);
/// }
/// ```
pub fn ifft<T>(data: &[T]) -> Vec<Complex>
where T: Into<Complex> + Copy {

    let length: isize = data.len() as isize;
    let norm: f64 = length as f64;

    // The pre-factor in the exponential terms
    let ipi_n: Complex = -Complex::i() * PI / length as f64;

    let mut res: Vec<Complex> = Vec::with_capacity(data.len());

    let mut an: Complex;
    let mut bkn: Complex;
    let mut bk_star: Complex;
    let mut sum_k: Complex;

    for k in 0..length {
        sum_k = Complex::default();
        bk_star = (k.pow(2) as f64 * ipi_n).exp().conjugate();

        for (n, val) in data.iter().enumerate() {
            an = (- (n.pow(2) as f64) * ipi_n).exp() * (*val).into();
            bkn = ((k - n as isize).pow(2) as f64 * ipi_n).exp();
            sum_k += an * bkn;
        }

        res.push(bk_star * sum_k / norm);
    }

    res
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
