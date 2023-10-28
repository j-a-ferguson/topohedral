//! This is the SPLline module
//! In it we compute:
//! - Non-Unform Rational Bsplines basis functions
//! 
//! 


use approx::ulps_eq;
use nalgebra as na;


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
const KNOT_ULPS: u32 = 32;

// --------------------------------------- Free Functions---------------------------------------- //

/// Tolerant less-thant for knots
fn knot_lt(u1: f64, u2: f64) -> bool 
{
    u1 < u2 || ulps_eq!(u1, u2, max_ulps = KNOT_ULPS) 
}

/// Tolerant greater-thant for knots
fn knot_gt(u1: f64, u2: f64) -> bool 
{
    u1 > u2 || ulps_eq!(u1, u2, max_ulps = KNOT_ULPS) 
}

/// Tolerant equal for knots
fn knot_eq(u1: f64, u2: f64) -> bool 
{
    ulps_eq!(u1, u2, max_ulps = KNOT_ULPS) 
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
        let end = (std::cmp::min(span_i, mi - 2 - pi)+1) as usize;
        let num_basis = end - start;
        (start, end, num_basis)
    }
    //..............................................................................................

    pub fn eval(&self, u: f64, p: usize, shape_funs: &mut [f64])
    {
        assert!(shape_funs.len() >= p+1, "Buffer too small to hold results");
        assert!(self.is_member(u), "u is outside of ");
        shape_funs.fill(0.0);
        shape_funs[0] = 1.0;

        let mut left = [0.0; PMAX];
        let mut right =[0.0; PMAX];

        let i = self.find_span(u, p);

        for j in 1..p+1
        {
            left[j-1] = u - self.knots[i + 1 - j];
            right[j-1] = self.knots[i + j] - u;
        }

        for j in 1..p+1
        {
            let mut saved = 0.0;
            for r in 0..j 
            {
                let temp = shape_funs[r] / (right[r] + left[j-r-1]);
                shape_funs[r] = saved + (right[r] * temp);
                saved = left[j - r - 1] * temp;
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
        assert!(shape_ders.len() >= (p+1)*(n+1), "Derivative buffer too small");

        shape_ders.fill(0.0);
        let mut shape_ders_mat = na::DMatrixViewMut::from_slice(shape_ders, n+1, p+1);

        let mut ndu_buf = [0.0; (PMAX+1) * (PMAX+1)];
        let mut ndu = na::DMatrixViewMut::from_slice(&mut ndu_buf, p+1, p+1);


        // ................... compute knot differences  ................... //
        let sp = self.find_span(u, p);
        let (mut saved, mut temp) = (0.0, 0.0);
        let mut left = [0.0; PMAX + 1];
        let mut right = [0.0; PMAX + 1];

        /* doc: explanation of ndu
        Upper triangle of ndu is the shape functions, so ndu[(0, 0)]
        contains N_{0,0}, ndu[(0, 1)] = N{0, 1}, ndu[(1, 1)] = N{1, 1} 
        and so on. 

        Lower triagle of ndu are the knot differences, so 
        ndu[(1,0)] = u_{i+1} - u_{i}
        ndu[(2,0)] = u_{i+1} - 
        */
        ndu[(0, 0)] = 1.0;

        for j in 1..p+1 
        {
            saved = 0.0;
            left[j] = u - self.knots[sp + 1 - j];
            right[j] = self.knots[sp + j] - u;

            for r in 0..j
            {
                ndu[(j, r)] = right[r + 1] + left[j - r];
                temp = ndu[(r, j - 1)] / ndu[(j, r)];

                ndu[(r, j)] = saved + (right[r+1] * temp);
                saved = left[j - r] * temp;
            }
            ndu[(j, j)] = saved;
        }



        // ................... compute derivatives ................... //

        // first column is the 0-order derivatives, so the shape funcitons themselves.
        shape_ders_mat.column_mut(0).copy_from(&ndu.column(p));
        let mut alpha_buf = [0.0; (PMAX+1)*(PMAX+1)];
        let mut alpha = na::DMatrixViewMut::from_slice(&mut alpha_buf, n+1, n+1);

        alpha[(0, 0)] = 1.0;
        for k in 1..n+1 
        {
            let u1 = self.knots[sp + p - k + 1];
            let u2 = self.knots[sp];
            alpha[(k, 0)] = alpha[(k-1, 0)] / (u1 - u2);

            for j in 1..k
            {
                let u3 = self.knots[sp + p + j - k + 1];
                let u4 = self.knots[sp + j];
                alpha[(k, j)] = (alpha[(k-1, j)] - alpha[(k-1, j-1)]) / (u3 - u4);
            }

            let u5 = self.knots[sp + p + 1];
            let u6 = self.knots[sp + k];
            alpha[(k, k)] = -alpha[(k-1, k-1)] / (u5 - u6);
        }

        for i in 0..p+1 // shape function index
        {
            for k in 1..n+1 // derivative order
            {
                for j in 0..k+1 // sum over shape funs
                {
                    shape_ders_mat[(i, k)] += alpha[(k, j)] * ndu[(j, p - k)];
                }

                let mut binom = (p - k + 1) as f64;
                for val in (p - k + 2) .. p+1 
                {
                    binom *= (val as f64);
                }
                shape_ders_mat[(i, k)] *= binom;
            }
        }


    }
}



// ------------------------------------------- Tests -------------------------------------------- //

#[cfg(test)]
mod tests {

    use crate::spl::KNOT_ULPS;

    use super::BsplineBasis;
    use super::PMAX;

    use approx::assert_relative_eq;
    use serde::Deserialize;
    use std::fs;
    use nalgebra as na;

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
                        assert_relative_eq!(val1, val2, max_relative = 1e-14);
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
    //..............................................................................................

    // macro_rules! eval_diff {
    //     ($test_name:ident, $knots:ident, $ders:ident, $order:expr) => {
    //         #[test]
    //         fn $test_name() {
    //             let test_data = TestData::new();
    //             let basis = BsplineBasis::new(test_data.$knots.values.clone());
    //             for (idx, u) in test_data.u.values.iter().enumerate()
    //             {
    //                 let ders_start = idx * ($order + 1);
    //                 let ders_end = (idx+1) * ($order + 1);
    //                 let ders1 = test_data.$ders.values[ders_start..ders_end].to_vec();


    //                 let mut ders_buf = [0.0; (PMAX+1) * (PMAX + 1)];
    //                 basis.eval_diff(*u, $order, $order, &mut ders_buf);
    //                 let ders2 = na::DMatrixViewMut::from_slice(&mut ders_buf, $order+1, $order+1);

                    
    //                 let (start, _end, num_basis) = basis.non_zero_basis(*u, $order);

    //                 // println!("........................................... {} {}", idx, u);
    //                 for j in 0..$order+1 
    //                 {
    //                     for k in 0..num_basis 
    //                     {
    //                         let val1 = ders1[j][k + start];
    //                         let val2 = ders2[(k, j)];
    //                         // println!("j = {} k = {} val1 = {} val2 = {}", j, k, val1, val2);
    //                         assert_relative_eq!(val1, val2, max_relative = 1e-14);
    //                     }
    //                 }
    //             }
    //         }
    //     };
    // }

    // eval_diff!(eval_diff0, knots_p0, ders_p0, 0);
    // eval_diff!(eval_diff1, knots_p1, ders_p1, 1);
    // eval_diff!(eval_diff2, knots_p2, ders_p2, 2);
    // eval_diff!(eval_diff3, knots_p3, ders_p3, 3);
    // eval_diff!(eval_diff4, knots_p4, ders_p4, 4);

}