/// # Boundary Enum
pub enum BoundaryD1<T>
where T: Fn(f64) -> f64 {
    FixedValue(f64),
    Function(T),
}

pub enum BoundaryDn<T>
where T: Fn(&[f64]) -> f64 {
    FixedValue(f64),
    Function(T),
}