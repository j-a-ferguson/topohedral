//! TODO add module docs here
//! 
//! 
//! 
//! 



use nalgebra as na;

use crate::spl;

use super::common::Curve;
use super::common::Vector;

// -------------------------------------------- Macros ------------------------------------------ //

// ------------------------------------------- Structs ------------------------------------------ //

pub struct Bcurve<const D: usize>
{
    p: usize, 
    basis: spl::BsplineBasis,
    points: Vec<Vector<D>>,
}

// ------------------------------------------- Enums -------------------------------------------- //


// ------------------------------------------- Types -------------------------------------------- //


// ----------------------------------------- Constants ------------------------------------------ //


// --------------------------------------- Free Functions---------------------------------------- //


// --------------------------------------- Implementations -------------------------------------- //

impl<const D: usize> Bcurve<D> {

    pub fn new<I1: IntoIterator<Item=f64>, 
               I2: IntoIterator<Item=Vector<D>>
               (p_in: usize, knots_in: I1, points_w_in: I2) -> Self 
    {
        Self { p: p_in, basis: (), points: () }
    }
}

impl<const D: usize> Curve for Bcurve<D>
{
    
    fn eval(&self, u: f64, point: &mut [f64])
    {
        let mut point_tmp = Vector::<D>::from_element(0.0);
        let (start, end, _nb) = self.basis.non_zero_basis(u, self.p);
        
        let mut basis_funs = [0.0; spl::PMAX];
        self.basis.eval(u, self.p, &mut basis_funs);

        for i in start .. end 
        {
            point_tmp += basis_funs[i - start] * self.points[i];
        }
        
        point.copy_from_slice(point_tmp.as_slice());
    }

    fn eval_diff(&self, u: f64, m: usize, diff: &mut [f64])
    {

    }

    fn eval_tangent(&self, u: f64, normalise: bool, tan: &mut [f64])
    {

    }

    fn eval_normal(&self, u: f64, normalise: bool, nor: &mut [f64])
    {


    }
    fn eval_binormal(&self, u: f64, normalise: bool, bin: &mut [f64])
    {


    }
    fn eval_curvature(&self, u: f64) -> f64
    {
        let out = 0.0;
        out

    }
    fn eval_torsion(&self, u: f64) -> f64
    {

        let out = 0.0;
        out

    }
    fn eval_arclen(&self, u1: f64, u2: f64) -> f64
    {


        let out = 0.0;
        out
    }
}

// ------------------------------------------- Tests -------------------------------------------- //

#[cfg(test)]
mod tests {

    use approx::assert_relative_eq;
    use serde::Deserialize;
    use std::fs;



    #[derive(Deserialize)]
    struct ParamData{
        description: String,
        values: Vec<f64>,
    }

    #[derive(Deserialize)]
    struct KnotData {
        description: String, 
        values: Vec<f64>,
    }

    #[derive(Deserialize)]
    struct WeightData {
        description: String, 
        values: Vec<f64>,
    }


    #[derive(Deserialize)]
    struct CpointData {
        description: String, 
        values: Vec<Vec<f64>>,
    }


    #[derive(Deserialize)]
    struct PointData {
        description: String, 
        values: Vec<Vec<f64>>,
    }

    #[derive(Deserialize)]
    struct DerData {
        description: String, 
        values: Vec<Vec<f64>>,
    }

    #[derive(Deserialize)]
    struct TangentData {
        description: String, 
        values: Vec<Vec<f64>>,
    }


    #[derive(Deserialize)]
    struct TestData {
        u: ParamData,
        knots_p1: KnotData,
        knots_p2: KnotData,
        knots_p3: KnotData,
        knots_p4: KnotData,
        weights_p1: WeightData,
        weights_p2: WeightData,
        weights_p3: WeightData,
        weights_p4: WeightData,
        cpoints_d2_p1: CpointData, 
        cpoints_d2_p2: CpointData, 
        cpoints_d2_p3: CpointData, 
        cpoints_d2_p4: CpointData, 
        cpoints_d3_p1: CpointData, 
        cpoints_d3_p2: CpointData, 
        cpoints_d3_p3: CpointData, 
        cpoints_d3_p4: CpointData, 
        points_d2_p1: PointData, 
        points_d2_p2: PointData, 
        points_d2_p3: PointData, 
        points_d2_p4: PointData, 
        points_d3_p1: PointData, 
        points_d3_p2: PointData, 
        points_d3_p3: PointData, 
        points_d3_p4: PointData, 
        ders_d2_p1: DerData, 
        ders_d2_p2: DerData, 
        ders_d2_p3: DerData, 
        ders_d2_p4: DerData, 
        ders_d3_p1: DerData, 
        ders_d3_p2: DerData, 
        ders_d3_p3: DerData, 
        ders_d3_p4: DerData, 
        tangent_d2_p1: TangentData, 
        tangent_d2_p2: TangentData, 
        tangent_d2_p3: TangentData, 
        tangent_d2_p4: TangentData, 
        tangent_d3_p1: TangentData, 
        tangent_d3_p2: TangentData, 
        tangent_d3_p3: TangentData, 
        tangent_d3_p4: TangentData, 
    }

    impl TestData {
        pub fn new() -> Self {
            let json_file = fs::read_to_string("assets/geo/bcurve-tests.json")
                                                        .expect("Unable to read file");
            serde_json::from_str(&json_file).expect("Could not deserialize")
        }
    }

    #[test]
    fn construction()
    {
        let test_data = TestData::new();
    }

}