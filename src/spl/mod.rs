//! This is the SPLline module
//! In it we compute:
//! - Non-Unform Rational Bsplines basis functions
//! 
//! 


use crate::utl;

use float_cmp::approx_eq;
use nalgebra as nl;


// ------------------------------------------- Structs ------------------------------------------ //

/// This struct represent a set of B-spline basis functions.
/// 
/// A B-spline basis is constructed from two things:
/// - A **knot vector** $\mathbf{u} = [u_{0}, ..., u_{m-1}]$ 
/// - A **polynomial degree** $p$
/// 
pub struct BsplineBasis
{
    knots: Vec<f64>,
}

// ------------------------------------------- Enums -------------------------------------------- //


// ----------------------------------------- Constants ------------------------------------------ //

/// This is the maximum allowable order of a bspline basis. It is an arbitrary number.
const PMAX: usize = 8;
/// This is the tolerance with which two knots are considered equal
const KNOT_ULPS: i64 = 256;

// --------------------------------------- Free Functions---------------------------------------- //

/// Tolerant less-thant for knots
fn knot_lt(u1: f64, u2: f64) -> bool 
{
    u1 < u2 || approx_eq!(f64, u1, u2, ulps = KNOT_ULPS) 
}

/// Tolerant greater-thant for knots
fn knot_gt(u1: f64, u2: f64) -> bool 
{
    u1 > u2 || approx_eq!(f64, u1, u2, ulps = KNOT_ULPS) 
}

/// Tolerant equal for knots
fn knot_eq(u1: f64, u2: f64) -> bool 
{
    approx_eq!(f64, u1, u2, ulps = KNOT_ULPS) 
}

fn upper_bound(arr: &[f64], value: f64) -> usize 
{
    let comp = |val: &f64| if knot_lt(*val, value)  {std::cmp::Ordering::Less} else {std::cmp::Ordering::Greater};
    arr.binary_search_by(comp).unwrap_or_else(|index| index)
}
// --------------------------------------- Implementations -------------------------------------- //
impl BsplineBasis
{
    pub fn new<I: IntoIterator<Item=f64>>(start: I) -> Self 
    {
        let knots: Vec<f64> = start.into_iter().collect();
        Self{knots}
    }
    //..............................................................................................

    pub fn is_member(&self, u: f64) -> bool
    {
        let umin = self.knots.first().unwrap();
        let umax = self.knots.last().unwrap();
        knot_gt(u, *umin) && knot_lt(u, *umax)
    }
    //..............................................................................................

    pub fn find_span(&self, u: f64, p: usize) -> usize 
    {
        let n = self.knots.len() - p - 1;
        let mut span: usize;

        if knot_eq(u, self.knots[n])
        {
            span = n-1;
        }
        else {
            let low = p;
            let idx = upper_bound(&self.knots[low..], u) + low;
            span = idx-1;
        }
        span
    }
    //..............................................................................................

    pub fn non_zero_basis(&self, u: f64, p: usize) -> (usize, usize, usize)
    {
        assert!(self.is_member(u));

        let mi = self.knots.len() as i32;
        let span_i = self.find_span(u, p) as i32;
        let pi = p as i32;

        let start = std::cmp::max(0, span_i - pi) as usize;
        let end = std::cmp::min(span_i, mi - 2 - pi) as usize;
        let num_basis = end - start;
        (start, end, num_basis)
    }
    //..............................................................................................

    pub fn eval(&self, u: f64, p: usize, shape_funs: &mut [f64])
    {
        shape_funs.fill(0.0);
        shape_funs[0] = 1.0;
        let i = self.find_span(u, p);
        let (mut saved, mut temp) = (0.0, 0.0);
        let (mut left, mut right) = ([0.0; PMAX], [0.0; PMAX]);    

        for j in 1..p+1
        {
            saved = 0.0;
            left[j] = u - self.knots[i + 1 - j];
            right[j] = self.knots[i + j] - u;

            for r in 0..j 
            {
                temp = shape_funs[r] / (right[r+1] + left[j-r]);
                shape_funs[r] = saved + (right[r+1] * temp);
                saved = left[j - r] * temp;
            }
            shape_funs[j] = saved;
        }
    }
    //..............................................................................................

    pub fn eval_diff(&self, u: f64, p: usize, n: usize, shape_ders: &mut [f64])
    {
        assert!(self.is_member(u), "u is not member of parameter range");
        assert!(p <= PMAX, "order exceeds max allowable");
        assert!(n <= p, "order of derivative exceeds order");


    }
}



// ------------------------------------------- Tests -------------------------------------------- //

#[cfg(test)]
mod tests {

    use crate::spl::KNOT_ULPS;

    use super::BsplineBasis;
    use super::PMAX;

    use float_cmp::assert_approx_eq;
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
    struct SpanData {
        description: String, 
        values: Vec<usize>,
    }

    #[derive(Deserialize)]
    struct BasisData {
        description: String, 
        values: Vec<Vec<f64>>,
    }


    #[derive(Deserialize)]
    struct TestData {
        u: ParamData, 
        knots_p0: KnotData, 
        knots_p1: KnotData, 
        knots_p2: KnotData, 
        knots_p3: KnotData, 
        knots_p4: KnotData, 
        span_p0: SpanData, 
        span_p1: SpanData, 
        span_p2: SpanData, 
        span_p3: SpanData, 
        span_p4: SpanData, 
        basis_p0: BasisData,
        basis_p1: BasisData,
        basis_p2: BasisData,
        basis_p3: BasisData,
        basis_p4: BasisData,
        ders_p0: BasisData,
        ders_p1: BasisData,
        ders_p2: BasisData,
        ders_p3: BasisData,
        ders_p4: BasisData,
    }

    impl TestData {
        pub fn new() -> Self
        {
            let json_file = fs::read_to_string("assets/spl/bsplinebasis.json").expect("Unable to read file");
            serde_json::from_str(&json_file).expect("Could not deserialize")
        }
    }


    #[test]
    fn construction()
    {
        let knots = vec![0.0, 0.1, 0.2, 0.3];
        let bspline_basis = BsplineBasis::new(knots.clone());
        assert_eq!(bspline_basis.knots, knots);

    }
    //..............................................................................................

    macro_rules! find_span {
        ($test_name:ident, $knots:ident, $span:ident, $order:expr) => {
            #[test]
            fn $test_name() {
                let test_data = TestData::new();
                let basis = BsplineBasis::new(test_data.$knots.values.clone());

                for (idx, u) in test_data.u.values.iter().enumerate()
                {
                    let span1 = test_data.$span.values[idx];
                    let span2 = basis.find_span(*u, $order);
                    if span1 != span2 
                    {
                        println!("idx = {} u = {} span1 = {} span2 = {}", idx, u, span1, span2);
                    }
                    assert_eq!(span1, span2);
                }
            }    
        };
    }
    find_span!(find_span0, knots_p0, span_p0, 0);
    find_span!(find_span1, knots_p1, span_p1, 1);
    find_span!(find_span2, knots_p2, span_p2, 2);
    find_span!(find_span3, knots_p3, span_p3, 3);
    find_span!(find_span4, knots_p4, span_p4, 4);
    //..............................................................................................

    macro_rules! eval {
        ($test_name:ident, $knots:ident, $basis:ident, $order:expr) => {
            #[test]
            fn $test_name() {

                let test_data = TestData::new();
                let basis = BsplineBasis::new(test_data.$knots.values.clone());

                for (idx, u) in test_data.u.values.iter().enumerate()
                {

                    let (start, end, _num_basis) = basis.non_zero_basis(*u, $order);
                    let basis_funs1 = test_data.$basis.values[idx].clone();

                    let mut basis_funs2 = [0.0; PMAX];
                    basis.eval(*u, $order, &mut basis_funs2);

                    for i in start..end
                    {
                        let val1 = basis_funs1[i];
                        let val2 = basis_funs2[i - start];
                        assert_approx_eq!(f64, val1, val2, ulps = KNOT_ULPS);
                    }
                    
                }
            }
        };
    }

    eval!(eval0, knots_p0, basis_p0, 0); 
    eval!(eval1, knots_p1, basis_p1, 1); 
    eval!(eval2, knots_p2, basis_p2, 2); 
    eval!(eval3, knots_p3, basis_p3, 3); 
    eval!(eval4, knots_p4, basis_p4, 4); 
    // eval!(eval0, knots_p0, basis_p0, 0); 

    // #[test]
    // fn eval() {

    //     let test_data = TestData::new();
    //     let basis = BsplineBasis::new(test_data.knots_p0.values.clone());

    //     for (idx, u) in test_data.u.values.iter().enumerate()
    //     {

    //         let (start, end, _num_basis) = basis.non_zero_basis(*u, 0);
    //         let basis_funs1 = test_data.basis_p0.values[idx].clone();

    //         let mut basis_funs2 = [0.0; PMAX];
    //         basis.eval(*u, 0, &mut basis_funs2);

    //         for i in start..end
    //         {
    //             let val1 = basis_funs1[i];
    //             let val2 = basis_funs2[i - start];
    //             assert_approx_eq!(f64, val1, val2, ulps = KNOT_ULPS);
    //         }
            
    //     }
    // }
}