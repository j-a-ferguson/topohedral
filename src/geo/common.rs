//! TODO add module docs here
//! 
//! 
//! 
//! 

use nalgebra as na;

// -------------------------------------------- Macros ------------------------------------------ //

// ------------------------------------------- Structs ------------------------------------------ //


// ------------------------------------------- Enums -------------------------------------------- //


// ------------------------------------------- Types -------------------------------------------- //

pub type Vector<const N: usize> = na::SVector<f64, N>;


// ------------------------------------------- Traits ------------------------------------------- //

/// This trait models the set of operations on a curve.
///  
pub trait Curve {

    fn eval(&self, u: f64, point: &mut [f64]);
    fn eval_diff(&self, u: f64, m: usize, diff: &mut [f64]);
    fn eval_diff_all(&self, u: f64, m: usize, diff: &mut[f64]);
    fn eval_tangent(&self, u: f64, normalise: bool, tan: &mut [f64]);
    fn eval_normal(&self, u: f64, normalise: bool, nor: &mut [f64]);
    fn eval_binormal(&self, u: f64, normalise: bool, bin: &mut [f64]);
    fn eval_curvature(&self, u: f64) -> f64;
    fn eval_torsion(&self, u: f64) -> f64;
    fn eval_arclen(&self, u1: f64, u2: f64) -> f64;
}

// ----------------------------------------- Constants ------------------------------------------ //


// --------------------------------------- Free Functions---------------------------------------- //

pub fn inv_homog<const N: usize>(point_w: &Vector<{N+1}>) -> Vector<{N}>
where 
    [(); N+1]:,
{
    let mut point = Vector::<{N}>::from_element(0.0);
    let w = point_w[N];
    for i in 0..N
    {
        point[i] = point_w[i] / w;
    }
    point
}

pub fn homog<const N: usize>(point: &Vector<N>, weight: f64) -> Vector<{N+1}>
{
    let mut point_w = Vector::<{N+1}>::from_element(0.0);
    for i in 0..N
    {
        point_w[i] = weight * point[i];
    }
    point_w[N] = weight;
    point_w
}

// --------------------------------------- Implementations -------------------------------------- //


// ------------------------------------------- Tests -------------------------------------------- //