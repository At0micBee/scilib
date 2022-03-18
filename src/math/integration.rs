//!
//! # Integration routines
//!

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

pub mod dx {

    fn chunk_quad(a: f64, b: f64, offset: f64) -> f64 {
        (a + b) / 2.0 * offset
    }

    fn chunk_trapez(a: f64, b: f64, offset: f64) -> f64 {
        if a.signum() != b.signum() {
            // Zero detected
            let a_abs = a.abs();
            let b_abs = b.abs();
            let thales = offset / (a_abs + b_abs);
            let first_tri_base = a_abs * thales;
            let second_tri_base = b_abs * thales;
            return (first_tri_base * a + second_tri_base * b) / 2.0;
        } else {
            let max;
            let min;
            if a < b {
                max = b;
                min = a;
            } else {
                max = a;
                min = b;
            }
            return (max - min) * offset / 2.0 + offset * min;
        }
    }

    pub fn compute_fn(
        function: impl Fn(f64) -> f64,
        lower_bound: f64,
        upper_bound: f64,
        div: usize,
        routine: fn(f64, f64, f64) -> f64,
    ) -> f64 {
        assert!(div >= 1);
        let mut area = 0.0;
        let offset = upper_bound - lower_bound;
        let mut a = function(lower_bound);
        let step = offset / div as f64;
        for i in 0..div {
            let b = function(step * (i as f64 + 1.0));
            area += routine(a, b, step);
            a = b;
        }
        area
    }

    pub fn compute_dt(data: &[(f64, f64)], routine: fn(f64, f64, f64) -> f64) -> f64 {
        let mut area = 0.0;
        for i in 0..data.len() - 1 {
            area += routine(data[i].1, data[i + 1].1, data[i + 1].0 - data[i].0);
        }
        area
    }

    #[inline]
    pub fn quad_fn(
        function: impl Fn(f64) -> f64,
        lower_bound: f64,
        upper_bound: f64,
        div: usize,
    ) -> f64 {
        compute_fn(function, lower_bound, upper_bound, div, chunk_quad)
    }

    #[inline]
    pub fn quad_dt(data: &[(f64, f64)]) -> f64 {
        compute_dt(data, chunk_quad)
    }

    #[inline]
    pub fn trapez_fn(
        function: impl Fn(f64) -> f64,
        lower_bound: f64,
        upper_bound: f64,
        div: usize,
    ) -> f64 {
        compute_fn(function, lower_bound, upper_bound, div, chunk_trapez)
    }

    #[inline]
    pub fn trapez_dt(data: &[(f64, f64)]) -> f64 {
        compute_dt(data, chunk_trapez)
    }
}

pub mod dxy {}

pub mod dxyz {}

pub mod dn {}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
