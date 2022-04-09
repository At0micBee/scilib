//! # Iterative methods to solve equation systems. 
//! 
//! Set of tools to and iterative methods to solve large space 
//! system of equations, using Krylov space based methods. 
//! This module contains : 
//!  - GMRES-Given  
//!  - Newton-Krylov 
//!  - Jacobian Free Newton-Krylov 

use std::{ops::{Add, Mul, Deref}, process::exit};


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


// /// # Matrix Vector multiplication 
// /// 
// /// Compute a matrix vector product
// /// 

// pub fn mat_vec_product<T>(mat : &Vec<Vec<T>> ,  vector : &Vec<T> )
// where T : Add + Mul {
    
//     assert_eq!(mat.len(),vector.len(),"Matrix and vector must have the same length")

//     mat.iter().zip(vector.iter()).map(||)
    
// }

/// # N-dimension Vector dot product
///  
/// Function computing the dot product of two 
/// vector in a N-dim space : 
/// $$ \sum_{i=1}^N u_i*v_i  $$
/// 
/// ```
/// # use scilib::math::krylov::dot_product;
/// 
/// let u = vec![1.0,4.0,1.0];
/// let v = vec![1.0,2.0,4.0];
/// 
/// let dot = dot_product(&u, &v);
/// 
/// assert_eq!(dot,13.0);
/// ```
/// 
pub fn dot_product<T>(vec1 : &Vec<T>, vec2 : &Vec<T>) -> T 
where T : Mul<T> + std::iter::Sum<<T as std::ops::Mul>::Output> + Copy{

    assert_eq!(vec1.len(),vec2.len()," Both vectors must have the same dimensions");

    vec1.iter().zip(vec2).map(|(v1,v2)| *v1 * *v2).sum()
}

/// # Jacobian Vector dot product estimate 
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

///  # Back Substitution Algorithm 
///  Function that solve the system 
///  $$ Ux = b $$  
///  where U is a upper triangular matrix (N,N), and b is a vector N. 
///  It use the following recursive relation 
///  $$ x_n = \frac{b_n}{U_{nn}} $$
///  $$ x_i =  \frac{b_i - \sum_{j=i+1}^n  U_{ij} x_j}{U_{ii}} $$
///   The algorithm can be found here <https://algowiki-project.org/en/Backward_substitution>
pub fn back_substitution(u: &Vec<Vec<f64>>,b:&Vec<f64>) -> Vec<f64> {

    // Size of the vector 
    let n = b.len();

    // Solution wich will be computed 
    let mut solution  = b.clone();

    // Initialize the loop 
    solution[n-1] = b[n-1]/u[n-1][n-1];

    // Do the sum 
    for j in n-1..0 {

        solution[j] = b[j]/u[j][j];

        for i in 0..j-1 {
            solution[i] = solution[i] - u[i][j] * solution[j];
        };
    };

    solution

}

/// # GMRES Given 
/// Given rotation GMRES version, it include the computation of the Arnoldi basis and 
/// the reduction of the Vk vector basis matrix. It return the du minimizing the 
/// residual (||func(u) - Jdu||₂) for a system of dimension N. Where func(u) is
/// the system we want to solve by JFNK, J is the jacobian of the system and du is  
/// the newton step we are looking for. For more information on the method :
/// Knoll et al 2009 and  Saad et al 1986.
pub fn gmres_given(
    func : fn(&Vec<f64>) -> Vec<f64>,   // function to minimize 
    u : Vec<f64>,                       // Point U where du is searched 
    du0 : Vec<f64>,                     // Initial guess  
    tol : f64,                          // Convergence tolerance 
    max_iter : u32,                     // Maximum number of iteration
) -> Vec<f64> {

    // Set up all needed matrix and vectors 

    let mut vk : Vec<Vec<f64>> = Vec::new() ;           // Arnoldi basis vector matrix 
    let mut vk_estimate : Vec<f64> = Vec::new() ;       // Intermediary krylov vector estimation 
    let mut hessian : Vec<Vec<f64>> = Vec::new();       // Hessenberg matrix, used to construct Vk    

    let mut norm_func : Vec<f64> = Vec::new() ;         // Evaluation of func(u) that will pass through the rotation
    let mut sn : Vec<f64> = Vec::new() ;                // Given rotation coefficients vectors
    let mut cn : Vec<f64> = Vec::new() ;                // Given rotation coefficients vectors
    let mut initial_fu : Vec<f64> = Vec::new() ;        // Evaluation of func(u)

    let mut residual : Vec<f64> = Vec::new() ;          // Residue vector 
    let mut residual_norm : f64  ;                      // Norm of the residue 

    // Evaluation of func(u)
    initial_fu = func(&u) ;
    
    // Compute initial residual and its norm 
    residual = jacobian_vec_estimate(&du0, &u, func); 
    residual_norm = l2_norm(&residual);

    // Compute the first krylov vector as Vk0 = residual/residual_norm
    vk.push(
        residual.iter().map(|r| r/residual_norm).collect()
    );

    norm_func.push(residual_norm); 


    // Beginning of the construction of the Krylov basis 
    while residual_norm > tol {

        // First estimation of the kth basis vector 
        vk_estimate = jacobian_vec_estimate(&vk[vk.len()], &u, func);


        exit(20)
    }; 

    initial_fu
    


}
