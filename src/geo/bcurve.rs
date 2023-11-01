//! TODO add module docs here
//! 
//! 
//! 
//! 



use nalgebra as na;

use crate::spl;

use super::common::Curve;

// -------------------------------------------- Macros ------------------------------------------ //

// ------------------------------------------- Structs ------------------------------------------ //

pub struct Bcurve<const D: usize>
{
    p: usize, 
    basis: spl::BsplineBasis,
    points: Vec<na::Vector3<f64>>,
}

// ------------------------------------------- Enums -------------------------------------------- //


// ------------------------------------------- Types -------------------------------------------- //


// ----------------------------------------- Constants ------------------------------------------ //


// --------------------------------------- Free Functions---------------------------------------- //


// --------------------------------------- Implementations -------------------------------------- //

impl<const D: usize> Curve for Bcurve<D>
{
    
    fn eval(u: f64, point: &mut [f64])
    {

    }

    fn eval_diff(u: f64, m: usize, diff: &mut [f64])
    {

    }

    fn eval_tangent(u: f64, normalise: bool, tan: &mut [f64])
    {

    }

    fn eval_normal(u: f64, normalise: bool, nor: &mut [f64])
    {


    }
    fn eval_binormal(u: f64, normalise: bool, bin: &mut [f64])
    {


    }
    fn eval_curvature(u: f64) -> f64
    {
        let out = 0.0;
        out

    }
    fn eval_torsion(u: f64) -> f64
    {

        let out = 0.0;
        out

    }
    fn eval_arclen(u1: f64, u2: f64) -> f64
    {


        let out = 0.0;
        out
    }
}

// ------------------------------------------- Tests -------------------------------------------- //