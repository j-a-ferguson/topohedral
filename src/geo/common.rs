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
    /// Evaluates a curve at the parameter value $u$
    /// # Arguments
    /// * `u` [I] - The curve parameter value
    /// * `point` [O] - slice of length $\geq D$
    fn eval(&self, u: f64, point: &mut [f64]);

    /// Evalutes the $m$'th derivative of the curve:
    /// $$
    ///     \mathbf{C}^{(m)}(u) = \frac{d^{m}C(u)}{du^{m}}
    /// $$
    /// # Arguments
    /// * `u` [I] - The curve parameter value
    /// * `ders` [O] - A slice of length $\geq D$. The i'th component 
    ///                is stored in `ders[i]`
    /// 
    fn eval_diff(&self, u: f64, m: usize, ders: &mut [f64]);

    /// Evalutes the $0$'th to the $m$'th derivative of the curve:
    /// $$
    ///     \left \\{ \mathbf{C}^{(0)}, ..., \mathbf{C}^{(m)}(u) \right \\}
    /// $$
    /// # Arguments
    /// * `u` [I] - The curve parameter value
    /// * `ders` [O] - A slice of length $\geq mD$. The i'th derivative  
    ///                is stored in `ders[i*D..(i+1)*D]`
    /// 
    fn eval_diff_all(&self, u: f64, m: usize, ders: &mut [f64]);
    fn eval_tangent(&self, u: f64, normalise: bool, tangent: &mut [f64]);
    fn eval_normal(&self, u: f64, normalise: bool, normal: &mut [f64]);
    fn eval_binormal(&self, u: f64, normalise: bool, binormal: &mut [f64]);
    fn eval_curvature(&self, u: f64) -> f64;
    fn eval_torsion(&self, u: f64) -> f64;
    fn eval_arclen(&self, u1: f64, u2: f64) -> f64;
}

// ----------------------------------------- Constants ------------------------------------------ //

// --------------------------------------- Free Functions---------------------------------------- //

pub fn inv_homog<const N: usize>(point_w: &Vector<{ N + 1 }>) -> Vector<{ N }>
where
    [(); N + 1]:,
{
    let mut point = Vector::<{ N }>::from_element(0.0);
    let w = point_w[N];
    for i in 0..N {
        point[i] = point_w[i] / w;
    }
    point
}

pub fn homog<const N: usize>(point: &Vector<N>, weight: f64) -> Vector<{ N + 1 }> {
    let mut point_w = Vector::<{ N + 1 }>::from_element(0.0);
    for i in 0..N {
        point_w[i] = weight * point[i];
    }
    point_w[N] = weight;
    point_w
}

// --------------------------------------- Implementations -------------------------------------- //

// ------------------------------------------- Tests -------------------------------------------- //
