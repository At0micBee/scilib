//!
//! # Basic math functions for series
//! 
//! This module provides access to many useful function that are not provided by the base Rust.
//! 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # Finds the maximum value in a slice
/// 
/// ## Definition
/// Looks for the largest number in the given set.
/// 
/// ## Inputs
/// - `val`: the slice with the numbers to compare
/// 
/// Returns the largest in the set.
/// 
/// ## Example
/// ```
/// # use scilib::math::series::max_slice;
/// let v: Vec<f64> = vec![0.0, 1.2, -0.1, 5.2, 0.254, 2.8];
/// let m: f64 = max_slice(&v);
/// assert_eq!(m, 5.2);
/// ```
pub fn max_slice<T>(val: &[T]) -> T 
where T: std::cmp::PartialOrd + Copy {

    let mut current_max: &T = &val[0];

    for v in val.iter() {
        if v > current_max {
            current_max = v;
        }
    }

    *current_max
}

/// # Finds the minimum value in a slice
/// 
/// ## Definition
/// Looks for the smallest number in the given set.
/// 
/// ## Inputs
/// - `val`: the slice with the numbers to compare
/// 
/// Returns the smallest in the set.
/// 
/// ## Example
/// ```
/// # use scilib::math::series::min_slice;
/// let v: Vec<f64> = vec![0.0, 1.2, -0.1, 5.2, 0.254, 2.8];
/// let m: f64 = min_slice(&v);
/// assert_eq!(m, -0.1);
/// ```
pub fn min_slice<T>(val: &[T]) -> T 
where T: std::cmp::PartialOrd + Copy {

    let mut current_min: &T = &val[0];

    for v in val.iter() {
        if v < current_min {
            current_min = v;
        }
    }

    *current_min
}

/// # Mean value of a series
/// 
/// ## Definition
/// We follow the mathematical definition of the mean:
/// $$
/// m = \frac{1}{n} \sum^{n}_{i = 1} x_i
/// $$
/// 
/// ## Inputs
/// - `val`: the slice of the series to compute
/// 
/// Returns the mean value of the series.
/// 
/// ## Example
/// ```
/// # use scilib::math::series::mean;
/// # use scilib::range;
/// let x: Vec<f64> = range::linear(0, 5, 6);
/// let m: f64 = mean(&x);
/// assert_eq!(m, 2.5);
/// ```
pub fn mean(val: &[f64]) -> f64 {

    val.iter().fold(0.0, |sum, v| sum + v) / val.len() as f64
}

/// # Variance of a series
/// 
/// ## Definition
/// The definition of the variance is:
/// $$
/// V = \frac{1}{n - 1} \sum^{n}_{i = 1} (x_i - m)^2
/// $$
/// Where $m$ is the mean of the series and $V$ the variance.
/// 
/// ## Inputs
/// - `val`: the slice of the series to compute
/// 
/// ## Example
/// ```
/// # use scilib::math::series::variance;
/// # use scilib::range;
/// let x: Vec<f64> = range::linear(0, 5, 6);
/// let v: f64 = variance(&x);
/// assert!((v - 3.5).abs() < 1e-10);
/// ```
pub fn variance(val: &[f64]) -> f64 {
    let mean: f64 = mean(val);
    val.iter().fold(0.0, |sum, v| sum + (v - mean).powi(2)) / (val.len() as f64 - 1.0)
}

/// # Standard deviation of a series
/// 
/// ## Definition
/// We follow the mathematical definition of the standard deviation:
/// $$
/// \sigma = \sqrt{ \frac{1}{n} \sum^{n}_{i = 1} (x_i - m)^2 } = \sqrt{V}
/// $$
/// Where $m$ is the mean of the series and $V$ the variance.
/// 
/// ## Inputs
/// - `val`: the slice of the series to compute
/// 
/// Returns the standard deviation value of the series.
/// 
/// ## Example
/// ```
/// # use scilib::math::series::std_dev;
/// # use scilib::range;
/// let x: Vec<f64> = range::linear(0, 5, 6);
/// let s: f64 = std_dev(&x);
/// assert!((s - 1.870828693386).abs() < 1e-10);
/// ```
pub fn std_dev(val: &[f64]) -> f64 {
    variance(val).sqrt()
}

/// # Skewness of a series
/// 
/// ## Definition
/// We follow the mathematical definition of the skewness:
/// $$
/// S = \frac{1}{n} \sum_{i-1}^{n} \left( \frac{x_i - m}{\sigma} \right)^3
/// $$
/// Where $m$ is the mean of the series and $\sigma$ the standard deviation.
/// 
/// ## Inputs
/// - `val`: the slice of the series to compute
/// 
/// ## Example
/// ```
/// # use scilib::math::series::skewness;
/// # use scilib::math::basic;
/// # use scilib::range;
/// let r: Vec<f64> = range::linear(-10, 10, 10000);
/// let g: Vec<f64> = r.iter().map(|x| basic::gaussian(1.0, 0.0, 1.7, *x)).collect();
/// let s: f64 = skewness(&g);
/// assert!((s - 1.348759).abs() <= 1e-3);
/// ```
pub fn skewness(val: &[f64]) -> f64 {
    let mean: f64 = mean(val);
    let sigma3: f64 = std_dev(val).powi(3);
    val.iter().fold(0.0, |sum, v| sum + (v - mean).powi(3)) / (val.len() as f64 * sigma3)
}

/// # Kurtosis of a series
/// 
/// ## Definition
/// We follow the mathematical definition of the kurtosis:
/// $$
/// S = \left[\frac{1}{n} \sum_{i-1}^{n} \left( \frac{x_i - m}{\sigma} \right)^4 \right] - 3
/// $$
/// Where $m$ is the mean of the series and $\sigma$ the standard deviation.
/// 
/// ## Inputs
/// - `val`: the slice of the series to compute
/// 
/// ## Example
/// ```
/// # use scilib::math::series::kurtosis;
/// # use scilib::math::basic;
/// # use scilib::range;
/// let r: Vec<f64> = range::linear(-10, 10, 10000);
/// let g: Vec<f64> = r.iter().map(|x| basic::gaussian(1.0, 0.0, 1.7, *x)).collect();
/// let s: f64 = kurtosis(&g);
/// assert!((s - 0.298867).abs() <= 1e-3);
/// ```
pub fn kurtosis(val: &[f64]) -> f64 {
    let mean: f64 = mean(val);
    let sigma4: f64 = std_dev(val).powi(4);
    val.iter().fold(0.0, |sum, v| sum + (v - mean).powi(4)) / (val.len() as f64 * sigma4) - 3.0
}

/// # Student's t value
/// 
/// ## Definition
/// Computes the difference of means between two series. The $t$ value is simply defined as:
/// $$
/// t = \frac{m_a - m_b}{s_D}
/// $$
/// where $m_a$ and $m_b$ are the means of series A and B, respectively. $s_D$ is the pooled invariance
/// of the two series, and is defined as:
/// $$
/// s_D=\sqrt{\frac{\sum_{a}(x_a-m_a)^2+\sum_{b}(x_b-m_b)^2}{N_a+N_b-2}\left(\frac{1}{N_a}+\frac{1}{N_b}\right)}
/// $$
/// where $x_a$ and $x_b$ are the points of their respective series, and $N_a$ and $N_b$ their number of points.
/// 
/// ## Inputs
/// - `val_a`: first series
/// - `val_b`: second series
/// 
/// ## Example
/// ```
/// # use scilib::math::series::student_t;
/// # use scilib::math::basic;
/// # use scilib::range;
/// let r: Vec<f64> = range::linear(-10, 10, 1000);
/// let g: Vec<f64> = r.iter().map(|x| basic::gaussian(1.0, 0.0, 1.7, *x)).collect();
/// let h: Vec<f64> = r.iter().map(|x| basic::gaussian(1.1, -0.2, 1.5, *x)).collect();
/// let t = student_t(&g, &h);
/// ```
pub fn student_t(val_a: &[f64], val_b: &[f64]) -> f64 {
    // First we compute the pooled variance
    let l_a: f64 = val_a.len() as f64;
    let l_b: f64 = val_b.len() as f64;
    let mean_a: f64 = mean(&val_a);
    let mean_b: f64 = mean(&val_b);
    let sum_a: f64 = val_a.iter().fold(0.0, |sum, v| sum + (v - mean_a).powi(2));
    let sum_b: f64 = val_b.iter().fold(0.0, |sum, v| sum + (v - mean_b).powi(2));
    let f1: f64 = (sum_a + sum_b) / (l_a + l_b - 2.0);
    let sd: f64 = (f1 * (1.0 / l_a + 1.0 / l_b)).sqrt();

    // Computing t
    (mean_a - mean_b) / sd
}

/// # Pearson r coefficient
/// 
/// ## Definition
/// The ![Pearson r coefficient](https://en.wikipedia.org/wiki/Pearson_correlation_coefficient)
/// is a correlation coefficient. Its use is widespread to check the correlation between two series
/// of data points. It is defined as:
/// $$
/// \rho_{X, Y} = \frac{\mathrm{cov}(X, Y)}{\sigma_X\sigma_Y}
/// = \frac{\sum_{i=0}^{n}(x_i - \bar x)(y_i - \bar y)}{\sqrt{\sum_{i=0}^{n} (x_i-\bar x)^2}\sqrt{\sum_{i=0}^{n} (y_i-\bar y)^2}}
/// $$
/// 
/// ## Inputs
/// - `sample_x`: the first series of values to check
/// - `sample_y`: the second series of values to check
/// 
/// Returns the Pearson r correlation coefficient between both series.
/// 
/// ## Example
pub fn pearson_r(sample_x: &[f64], sample_y: &[f64]) -> f64 {
    
    let mean_x: f64 = mean(sample_x);   // Computing mean for x
    let mean_y: f64 = mean(sample_y);   // Computing mean for y

    let mut temp_x: f64;                // Creating temporary value for x
    let mut temp_y: f64;                // Creating temporary value for y
    let mut t: f64 = 0.0;               // Top part of Pearson
    let mut b_x: f64 = 0.0;             // First div of Pearson
    let mut b_y: f64 = 0.0;             // Second div of Pearson

    for (x, y) in sample_x.iter().zip(sample_y) {
        temp_x = x - mean_x;
        temp_y = y - mean_y;
        t += temp_x * temp_y;
        b_x += temp_x.powi(2);
        b_y += temp_y.powi(2);
    }

    t / (b_x * b_y).sqrt()
}

/// # Min-Max scaling of a series
/// 
/// ## Definition
/// Min-Max scaling compresses all the data points passed in a series between two arbitrary values a and b.
/// $$
/// x_{s} = \frac{a + (x - min(x))(b - a)}{max(x) - min(x)}
/// $$
/// 
/// ## Inputs
/// - `val`: the series to scale
/// - `a`: the minimum to scale to
/// - `b`: the maximum to scale to
/// 
/// Returns the new series between a and b.
/// 
/// ## Example
/// ```
/// # use scilib::range;
/// # use scilib::math::series::scale_min_max;
/// let x: Vec<f64> = range::linear(1, 6, 7);
/// let n: Vec<f64> = scale_min_max(&x, 2.0, -1.0);
/// assert_eq!(n[0], 2.0);
/// assert_eq!(n[3], 0.5);
/// assert_eq!(n[6], -1.0);
/// ```
pub fn scale_min_max(val: &[f64], a: f64, b: f64) -> Vec<f64> {

    let max_val: f64 = max_slice(&val);     // We find the max of the slice
    let min_val: f64 = min_slice(&val);     // We find the min of the slice
    let ba: f64 = b - a;                    // We compute the difference for the top
    let div: f64 = max_val - min_val;       // We compute the difference for the divisor

    val.iter().map(|x| {                    // We go through each value to scale
        a + (x - min_val) * ba / div        // Scaling
    }).collect()                            // Returning the right type of vector
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
