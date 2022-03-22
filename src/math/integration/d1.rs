fn chunk_quad(a: f64, offset: f64) -> f64 {
    a * offset
}

fn chunk_trapz(a: f64, b: f64, offset: f64) -> f64 {
    (a + b) / 2.0 * offset
}

fn chunk_simp(
    p0: (f64 /*a0*/, f64 /*b0*/),
    p1: (f64 /*a1*/, f64 /*b1*/),
    p2: (f64 /*a2*/, f64 /*b2*/),
) -> f64 {
    let coeff0 = p0.1 / ((p0.0 - p1.0) * (p0.0 - p2.0));
    let coeff1 = p1.1 / ((p1.0 - p0.0) * (p1.0 - p2.0));
    let coeff2 = p2.1 / ((p2.0 - p0.0) * (p2.0 - p1.0));
    let a1_a2 = p1.0 + p2.0;
    let a0_a2 = p0.0 + p2.0;
    let a0_a1 = p0.0 + p1.0;
    let a1a2 = p1.0 * p2.0;
    let a0a2 = p1.0 * p2.0;
    let a0a1 = p0.0 * p1.0;
    let l0 = coeff0 * (p0.0.powi(3) / 3.0 - p0.0.powi(2) / 2.0 * a1_a2 + p0.0 * a1a2)
        + coeff1 * (p0.0.powi(3) / 3.0 - p0.0.powi(2) / 2.0 * a0_a2 + p0.0 * a0a2)
        + coeff2 * (p0.0.powi(3) / 3.0 - p0.0.powi(2) / 2.0 * a0_a1 + p0.0 * a0a1);
    let l2 = coeff0 * (p2.0.powi(3) / 3.0 - p2.0.powi(2) / 2.0 * a1_a2 + p2.0 * a1a2)
        + coeff1 * (p2.0.powi(3) / 3.0 - p2.0.powi(2) / 2.0 * a0_a2 + p2.0 * a0a2)
        + coeff2 * (p2.0.powi(3) / 3.0 - p2.0.powi(2) / 2.0 * a0_a1 + p2.0 * a0a1);
    l2 - l0
}

pub fn quad_fn(
    function: impl Fn(f64) -> f64,
    lower_bound: f64,
    upper_bound: f64,
    div: usize,
) -> f64 {
    assert!(div >= 1);
    let mut area = 0.0;
    let step = (upper_bound - lower_bound) / div as f64;
    for i in 0..div {
        area += chunk_quad(function(step * (i as f64 + 1.0) + lower_bound), step);
    }
    area
}

pub fn trapz_fn(
    function: impl Fn(f64) -> f64,
    lower_bound: f64,
    upper_bound: f64,
    div: usize,
) -> f64 {
    assert!(div >= 1);
    let mut area = 0.0;
    let step = (upper_bound - lower_bound) / div as f64;
    let mut a = function(lower_bound);
    for i in 0..div {
        let b = function(step * (i as f64 + 1.0) + lower_bound);
        area += chunk_trapz(a, b, step);
        a = b;
    }
    area
}

pub fn simp_fn(
    function: impl Fn(f64) -> f64,
    lower_bound: f64,
    upper_bound: f64,
    div: usize,
) -> f64 {
    assert!(div >= 2 && div % 2 == 0);
    let mut area = 0.0;
    let step = (upper_bound - lower_bound) / div as f64;
    let mut p0 = (lower_bound, function(lower_bound));
    let mut curr_pos = lower_bound + step;
    let mut p1 = (curr_pos, function(curr_pos));
    for i in (1..div).step_by(2) {
        curr_pos = lower_bound + (i + 2) as f64 * step;
        let p2 = (curr_pos, function(curr_pos));
        area += chunk_simp(p0, p1, p2);
        p0 = p1;
        p1 = p2;
    }
    area
}

pub fn quad_dt_tup(data: &[(f64, f64)]) -> f64 {
    let mut area = 0.0;
    for i in 0..data.len() - 1 {
        area += chunk_quad(data[i].1, data[i + 1].0 - data[i].0);
    }
    area
}

pub fn trapz_dt_tup(data: &[(f64, f64)]) -> f64 {
    let mut area = 0.0;
    for i in 0..data.len() - 1 {
        area += chunk_trapz(data[i].1, data[i + 1].1, data[i + 1].0 - data[i].0);
    }
    area
}

pub fn simp_dt_tup(data: &[(f64, f64)]) -> f64 {
    assert!(data.len() >= 3 && data.len() % 2 == 1);
    let mut area = 0.0;
    for i in (0..data.len() - 2).step_by(2) {
        area += chunk_simp(data[i], data[i + 1], data[i + 2]);
    }
    area
}

pub fn quad_dt_sep(x: &[f64], y: &[f64]) -> f64 {
    assert!(x.len() == y.len());
    let mut area = 0.0;
    for i in 0..x.len() - 1 {
        area += chunk_quad(y[i], x[i + 1] - x[i]);
    }
    area
}

pub fn trapz_dt_sep(x: &[f64], y: &[f64]) -> f64 {
    assert!(x.len() == y.len());
    let mut area = 0.0;
    for i in 0..x.len() - 1 {
        area += chunk_trapz(y[i], y[i + 1], x[i + 1] - x[i]);
    }
    area
}

pub fn simp_dt_sep(x: &[f64], y: &[f64]) -> f64 {
    assert!(x.len() == y.len() && x.len() >= 3 && x.len() % 2 == 1);
    let mut area = 0.0;
    for i in (0..x.len() - 2).step_by(2) {
        area += chunk_simp((x[i], y[i]), (x[i + 1], y[i + 1]), (x[i + 2], y[i + 2]));
    }
    area
}
