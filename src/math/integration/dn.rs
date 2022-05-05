////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// # (dn) Rectangle rule
/// Integrate n-dimensional function through the rectangle rule:
pub fn rect(
    function: impl Fn(&[f64]) -> f64,
    upper_boundaries: &[impl Fn(&[f64]) -> f64],
    lower_boundaries: &[impl Fn(&[f64]) -> f64],
    divs: &[usize]
) -> f64 {

}

/// # (dn) Trapezoidal rule
/// Integrate n-dimensional function through the trapezoidal rule:
pub fn trapz(
    function: impl Fn(&[f64]) -> f64,
    upper_boundaries: &[impl Fn(&[f64]) -> f64],
    lower_boundaries: &[impl Fn(&[f64]) -> f64],
    divs: &[usize]
) -> f64 {

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
