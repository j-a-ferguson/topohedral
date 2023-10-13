
use nalgebra as na;

use crate::spl;

pub struct Bcurve
{
    points: Vec<na::Vector3<f64>>,
    basis: spl::BsplineBasis
}