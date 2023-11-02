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

pub type Vector<const D: usize> = na::SVector<f64, D>;

// ------------------------------------------- Traits ------------------------------------------- //

/// This trait models the set of operations on a curve.
///  
pub trait Curve {

    fn eval(&self, u: f64, point: &mut [f64]);
    fn eval_diff(&self, u: f64, m: usize, diff: &mut [f64]);
    fn eval_tangent(&self, u: f64, normalise: bool, tan: &mut [f64]);
    fn eval_normal(&self, u: f64, normalise: bool, nor: &mut [f64]);
    fn eval_binormal(&self, u: f64, normalise: bool, bin: &mut [f64]);
    fn eval_curvature(&self, u: f64) -> f64;
    fn eval_torsion(&self, u: f64) -> f64;
    fn eval_arclen(&self, u1: f64, u2: f64) -> f64;
}

// ----------------------------------------- Constants ------------------------------------------ //


// --------------------------------------- Free Functions---------------------------------------- //

pub fn inv_homog<const D: usize>(point_w: &Vector<D>) -> Vector<{D-1}>
{
    let mut point = Vector::<{D-1}>::from_element(0.0);
    let w = point_w[D-1];
    for i in 0..D-1
    {
        point[i] = point_w[i] / w;
    }
    point
}

pub fn homog<const D: usize>(point: &Vector<D>, weight: f64) -> Vector<{D+1}>
{
    let mut point_w = Vector::<{D+1}>::from_element(0.0);
    for i in 0..D
    {
        point_w[i] = weight * point[i];
    }
    point_w[D] = weight;
    point_w
}

// --------------------------------------- Implementations -------------------------------------- //


// ------------------------------------------- Tests -------------------------------------------- //