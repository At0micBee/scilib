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
pub fn convolve(a: &Vec<f64>, b: &Vec<f64>) -> Vec<f64> {

    let mut sum: f64;
    let l_a: usize = a.len();
    let l_b: usize = b.len();
    let mut res: Vec<f64> = Vec::with_capacity(a.len());

    // Box b starts to slide over a
    for n in 0..(l_b - 1) {
        sum = 0.0;
        for m in 0..=n {
            sum += a[n-m] * b[m];

        }
        res.push(sum);
    }

    // Box b covers a completely
    for n in (l_b - 1)..=(l_a - 1) {
        sum = 0.0;
        for m in 0..l_b {
            sum += a[n-m] * b[m];
        }
        res.push(sum);
    }

    // Box b exits the slide over a
    for n in l_a..(l_a + l_b - 1) {
        sum = 0.0;
        for m in (n - l_a + 1)..l_b {
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
pub fn fft<T>(data: &Vec<T>) -> Vec<Complex>
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
        sum_k = 0.into();
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
pub fn ifft<T>(data: &Vec<T>) -> Vec<Complex>
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
        sum_k = 0.into();
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
