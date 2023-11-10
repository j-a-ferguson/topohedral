//! TODO add module docs here
//! 
//! 
//! 
//! 




use crate::spl;
use crate::geo::common::{Vector, Curve, inv_homog, homog};
use crate::utl::NDArray;

// -------------------------------------------- Macros ------------------------------------------ //

// ------------------------------------------- Structs ------------------------------------------ //

/// This struct encapsulates the data for a B-curve
/// 
pub struct Bcurve<const D: usize>
where 
    [(); D+1]:,
{
    p: usize, 
    knots: Vec<f64>,
    cpoints_w: Vec<Vector<{D+1}>>
}

// ------------------------------------------- Enums -------------------------------------------- //


// ------------------------------------------- Types -------------------------------------------- //


// ----------------------------------------- Constants ------------------------------------------ //


// --------------------------------------- Free Functions---------------------------------------- //

// fn binom_coeff(p: usize, k: usize) -> f64
// {
//     assert!(p >= k, "p must be greater than k");
//     let mut out = (p - k + 1);
//     for i in p-k+2 .. p+1
//     {
//         out *= i;
//     }
//     out as f64
// }

fn binom_coeff(n: usize, binom:  &mut [f64])
{
    debug_assert!(binom.len() >= (n+1)*(n+1) ) ;

    binom.fill(0.0);
    let mut binom_arr = NDArray::<'_, f64, 2>::new(binom, &[n+1, n+1]);

    for i in 0..n+1
    {
        binom_arr[&[i, i]] = 1.0;
        binom_arr[&[i, 0]] = 1.0;
    }

    for n2 in 2..n+1
    {
        for k2 in 1..n2
        {
            binom_arr[&[n2, k2]] = binom_arr[&[n2-1, k2 - 1]] + binom_arr[&[n2-1, k2]];
        }
    }
} 

// --------------------------------------- Implementations -------------------------------------- //

impl<const D: usize> Bcurve<D> 
where 
    [(); D+1]:,
{

    pub fn new(p: usize, knots: &[f64], weights: &[f64], cpoints: &[Vector<D>]) -> Self
    {
        debug_assert!(p <= spl::PMAX, "Order too large");
        debug_assert!(knots.is_sorted(), "knots not sorted");
        debug_assert!(weights.iter().all(|&x| x >= 0.0));
        debug_assert!(weights.len() == cpoints.len());
        debug_assert!(knots.len() == cpoints.len() + p + 1);

        let mut points_w = vec![Vector::<{D+1}>::zeros(); cpoints.len()];

        for i in 0..cpoints.len()
        {
            points_w[i] = homog(&cpoints[i], weights[i]);
        }

        Self {
            p: p, 
            knots: knots.to_vec(), 
            cpoints_w: points_w,
        }
    }
}

impl<const D: usize> Curve for Bcurve<D>
where 
    [(); D+1]:,
{
    
    fn eval(&self, u: f64, point: &mut [f64])
    {
        debug_assert!(spl::is_member(&self.knots, u));
        debug_assert!(point.len() >= D);

        let mut pointw_tmp = Vector::<{D+1}>::from_element(0.0);
        let (start, end, _nb) = spl::non_zero_basis(&self.knots, u, self.p);
        
        let mut basis_funs = [0.0; spl::PMAX];
        spl::eval(&self.knots, u, self.p, &mut basis_funs);

        for i in start .. end 
        {
            pointw_tmp += basis_funs[i - start] * self.cpoints_w[i];
        }

        let point_tmp = inv_homog(&pointw_tmp);
        point.copy_from_slice(point_tmp.as_slice());
    }

    fn eval_diff(&self, u: f64, m: usize, diff: &mut [f64])
    {
        debug_assert!(spl::is_member(&self.knots, u));
        debug_assert!(diff.len() >= D);

        let mut dersw = [Vector::<{D+1}>::zeros(); spl::PMAX];
        let (start, end, num_basis) = spl::non_zero_basis(&self.knots, u, self.p);

    }

    fn eval_diff_all(&self, u: f64, k: usize, ders: &mut[f64])
    {
        debug_assert!(spl::is_member(&self.knots, u));
        debug_assert!(ders.len() >= D);
        debug_assert!(k <= self.p);

        const DIM_MAX: usize = spl::PMAX+1;
        let dim = k+1;
        let mut dersw = [Vector::<{D+1}>::zeros(); DIM_MAX];
        let (start, _, num_basis) = spl::non_zero_basis(&self.knots, u, self.p);

        let mut basis_ders = [0.0; DIM_MAX * DIM_MAX];
        spl::eval_diff_all(&self.knots, u, self.p, k, &mut basis_ders);
        let basis_ders_arr  = NDArray::<'_, f64, 2>::new(&mut basis_ders, &[dim+1, dim+1]);

        for m in 0..k+1 // loop over derivatives
        {
            for j in 0..num_basis
            {
                let nj = basis_ders_arr[&[j, m]];
                let pwj = self.cpoints_w[start + j];
                dersw[m] += nj * pwj;
            }
        }

        let mut binom = [0.0; DIM_MAX * DIM_MAX];
        binom_coeff(k+1, &mut binom);
        let binom_arr = NDArray::<'_, f64, 2>::new(&mut binom, &[dim, dim]);

        for m in 0..k+1
        {
            let mut v = dersw[m].clone();
            for j in 1..m+1
            {
                v -= binom_arr[&[m, j]] * dersw[m - j];
            }
        }
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

    use approx::{assert_relative_eq, ulps_eq};
    use serde::Deserialize;
    use std::fs;

    use crate::geo::common::Curve;
    use crate::utl::NDArray;

    use super::Bcurve;
    use super::Vector;
    use super::binom_coeff;

    
    fn convert<const D: usize>(data: &Vec::<Vec<f64>>) -> Vec<Vector<D>>
    {
        for val in data 
        {
            assert!(val.len() == D);
        }

        // let mut out = Vec::<Vector<D>>::with_capacity(data.len());
        let mut out = vec![Vector::<D>::zeros(); data.len()];

        for (idx, val) in data.iter().enumerate()
        {
            out[idx].copy_from_slice(val);
        }
        out
    }



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
    fn binomcoeff()
    {
        let mut binom = [0.0; 36];
        binom_coeff(5, &mut binom);
        let binom_arr = NDArray::<'_, f64, 2>::new(&mut binom, &[6,6]);
        assert!(ulps_eq!(binom_arr[&[5, 0]], 1.0, max_ulps = 4));
        assert!(ulps_eq!(binom_arr[&[5, 1]], 5.0, max_ulps = 4));
        assert!(ulps_eq!(binom_arr[&[5, 2]], 10.0, max_ulps = 4));
        assert!(ulps_eq!(binom_arr[&[5, 3]], 10.0, max_ulps = 4));
        assert!(ulps_eq!(binom_arr[&[5, 4]], 5.0, max_ulps = 4));
        assert!(ulps_eq!(binom_arr[&[5, 5]], 1.0, max_ulps = 4));
    }

    #[test]
    fn construction()
    {
        let test_data = TestData::new();
        let p = 3;
        let knots = test_data.knots_p3.values;
        let weights = test_data.weights_p3.values;
        let cpoints: Vec<Vector<2>> = convert(&test_data.cpoints_d2_p3.values);
        let bcurve = Bcurve::<2>::new(p, &knots, &weights, &cpoints);

    }
    //..............................................................................................

    macro_rules! eval {
        ($test_name: ident, $knots: ident, $weights: ident, $cpoints:ident, $points: ident, $dim: expr, $order:expr) => {
            #[test]
            fn $test_name() {
                let d = $dim;
                let p = $order;
                let test_data = TestData::new();
                let knots = test_data.$knots.values;
                let weights = test_data.$weights.values;
                let cpoints: Vec<Vector<$dim>> = convert(&test_data.$cpoints.values);
                let points = test_data.$points.values;
                let bcurve = Bcurve::<$dim>::new(p, &knots, &weights, &cpoints);

                for (idx, u) in test_data.u.values.iter().enumerate()
                {
                    let point1 = points[idx].clone();
                    let mut point2 = Vector::<$dim>::zeros();
                    bcurve.eval(*u, point2.as_mut_slice()); 
                    for i in 0..d
                    {
                        assert_relative_eq!(point1[i], point2[i], max_relative = 1e-12);
                    }
                }
            }     
        };
    }
    eval!(eval_d2_p1, knots_p1, weights_p1, cpoints_d2_p1, points_d2_p1, 2, 1);
    eval!(eval_d2_p2, knots_p2, weights_p2, cpoints_d2_p2, points_d2_p2, 2, 2);
    eval!(eval_d2_p3, knots_p3, weights_p3, cpoints_d2_p3, points_d2_p3, 2, 3);
    eval!(eval_d2_p4, knots_p4, weights_p4, cpoints_d2_p4, points_d2_p4, 2, 4);
    eval!(eval_d3_p1, knots_p1, weights_p1, cpoints_d3_p1, points_d3_p1, 3, 1);
    eval!(eval_d3_p2, knots_p2, weights_p2, cpoints_d3_p2, points_d3_p2, 3, 2);
    eval!(eval_d3_p3, knots_p3, weights_p3, cpoints_d3_p3, points_d3_p3, 3, 3);
    eval!(eval_d3_p4, knots_p4, weights_p4, cpoints_d3_p4, points_d3_p4, 3, 4);
    //..............................................................................................

    macro_rules! eval_diff {
        ($test_name: ident, $knots: ident, $weights: ident, $cpoints:ident, $ders: ident, $dim: expr, $order:expr) => {
            #[test]
            fn $test_name() {
                let d = $dim;
                let p = $order;
                let test_data = TestData::new();
                let knots = test_data.$knots.values;
                let weights = test_data.$weights.values;
                let cpoints: Vec<Vector<$dim>> = convert(&test_data.$cpoints.values);
                let ders = test_data.$ders.values;
                let bcurve = Bcurve::<$dim>::new(p, &knots, &weights, &cpoints);

                for (idx, u) in test_data.u.values.iter().enumerate()
                {
                    let ders_all_1 = ders[idx].clone();
                    for k in 0..p+1
                    {
                        let ders1 = &ders_all_1[(k * d) .. ((k + 1) * d)];
                        let mut ders2 = Vector::<$dim>::zeros();
                        bcurve.eval_diff(*u, k, ders2.as_mut_slice());
                        for i in 0..d
                        {
                            assert_relative_eq!(ders1[i], ders2[i], max_relative = 1e-12);
                        }
                    }
                }
            }     
        };
    }
    eval_diff!(eval_diff_d2_p1, knots_p1, weights_p1, cpoints_d2_p1, ders_d2_p1, 2, 1);
    // eval!(eval_d2_p2, knots_p2, weights_p2, cpoints_d2_p2, points_d2_p2, 2, 2);
    // eval!(eval_d2_p3, knots_p3, weights_p3, cpoints_d2_p3, points_d2_p3, 2, 3);
    // eval!(eval_d2_p4, knots_p4, weights_p4, cpoints_d2_p4, points_d2_p4, 2, 4);
    // eval!(eval_d3_p1, knots_p1, weights_p1, cpoints_d3_p1, points_d3_p1, 3, 1);
    // eval!(eval_d3_p2, knots_p2, weights_p2, cpoints_d3_p2, points_d3_p2, 3, 2);
    // eval!(eval_d3_p3, knots_p3, weights_p3, cpoints_d3_p3, points_d3_p3, 3, 3);
    // eval!(eval_d3_p4, knots_p4, weights_p4, cpoints_d3_p4, points_d3_p4, 3, 4);

}