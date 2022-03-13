//! # Iterative methods to solve equation systems. 
//! 
//! Set of tools to and iterative methods to solve large space 
//! system of equations, using Krylov space based methods. 
//! This module contains : 
//!  - GMRES-Given  
//!  - Newton-Krylov 
//!  - Jacobian Free Newton-Krylov 

// Test function for all run test
pub fn test_func(v : &Vec<f64>) -> Vec<f64> {

    let res_1 = v[0] - 3.0* v[1] + v[2] - 2.0;
    let res_2 = 3.0* v[0] - 4.0 * v[1] + v[2]; 
    let res_3 = 2.0* v[1] - v[2] - 1.0;

    let res_vec = vec![res_1,res_2,res_3];
    
    res_vec
}

/// # L2 Norm 
/// 
///  Compute L2 norm of a N dim vector, computed as 
/// $$ |V| = \sqrt{ \sum_{i=0}^N v_i^2 } $$ 
/// ```
/// # use::scilib::math::krylov::l2_norm;
/// let v = vec![2.0_f64.sqrt(),2.0f64.sqrt()];
/// 
/// assert!(l2_norm(&v) == 2.0)
/// ```
pub fn l2_norm(u : &Vec<f64>) -> f64 {

    // compute the norm 
    let norm : f64 = u.iter(). // Create iterator 
    map(|x : &f64| x*x) // x*x 
    .sum::<f64>() // sum x_i^2
    .sqrt(); // sqrt(sum x_i^2)

    norm
}

/// Estimation of the dot product J.v of the equation system
/// J.v is estimated at point U as:
///  $$ \frac{F(U+\epsilon V) - F(U)}{\epsilon}$$
/// 
/// ```
/// # use scilib::math::krylov::{jacobian_vec_estimate,test_func};
/// 
/// let u = vec![1.0,1.0,1.0];
/// let v = vec![1e-1,1e-1,1e-1];
///
/// let jv = jacobian_vec_estimate(&u, &v, test_func);
///
/// assert!((jv[0] + 1.0).abs() < 1e-4 && jv[1].abs() < 1e-4 && (jv[2] - 1.0).abs()<1e-4);
/// 
/// ```
pub fn jacobian_vec_estimate(
    v : &Vec<f64>, 
    u : &Vec<f64>, 
    func : fn(&Vec<f64>) -> Vec<f64> ) -> Vec<f64> {
    
    let norm_u : f64 = l2_norm(u); // |U|
    let norm_v : f64 = l2_norm(v); // |V|
    
    let mut jac_vec : Vec<f64> = vec![0.0;u.len()];
    
    // If the perturbation vector have a null norm then the dot product is a null vector
    // It then stay full of 0.0
    if norm_v != 0.0 {

        // Adapted epsilon to estimate the J.v dot product 
        let eps : f64 = ((1.0 + norm_u) * f64::EPSILON ).sqrt() / norm_v ; 

        // F(U)
        let func_u : Vec<f64> = func(&u);
        
        // W = U + eps*V
        let w : Vec<f64> = u.iter().zip(v.iter())
        .map(|(x,y)| x + eps*y ).collect();

        // F(W)
        let func_eps_v = func(&w);

        // J.v = (F(W) - F(U)) / eps
        jac_vec = func_eps_v.iter().zip(func_u.iter())
        .map(|(x,y)| (x - y)/eps).collect();
    }

    jac_vec
}


