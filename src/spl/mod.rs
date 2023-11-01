//! This is the SPLline module
//! In it we compute:
//! - Non-Unform Rational Bsplines basis functions
//! 
//! 


use approx::ulps_eq;
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

fn binom_coeff(p: usize, k: usize) -> f64
{
    assert!(p >= k, "p must be greater than k");
    let mut out = (p - k + 1);
    for i in p-k+2 .. p+1
    {
        out *= i;
    }
    out as f64
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
            left[j-1] = u - self.knots[i - p + j];
            right[j-1] = self.knots[i + j] - u;
        }

        for j in 1..p+1
        {
            let mut saved = 0.0;
            for r in 0..j 
            {
                let ri = right[r];
                let le = left[p - j + r];
                let temp = shape_funs[r] / (ri + le);
                shape_funs[r] = saved + (ri * temp);
                saved = le * temp;
            }
            shape_funs[j] = saved;
        }
    }
    //..............................................................................................



    pub fn eval_diff(&self, u: f64, p: usize, k: usize, shape_ders: &mut [f64])
    {
        assert!(self.is_member(u), "u is not member of parameter range");
        assert!(p <= PMAX, "order exceeds max allowable");
        assert!(k <= p, "order of derivative exceeds order");
        assert!(shape_ders.len() >= (p+1), "Derivative buffer too small");

        let i = self.find_span(u, p);
        let p_k = p - k;

        /* doc: compute shape functions
        This code-block initialises the first p_k+1 elements of 
        shape_ders to the shape functions: 

            N_{i-p_k, p_k}^{(0)} ... N_{i-p_k, p_k}^{(0)}
        
        This is the first phase of the cox-deboor recursion.
         */
        {
            shape_ders.fill(0.0);
            shape_ders[0] = 1.0;
            let mut left = [0.0; PMAX];
            let mut right =[0.0; PMAX];
            for j in 1..p_k+1
            {
                left[j-1] = u - self.knots[i - p_k + j];
                right[j-1] = self.knots[i + j] - u;
            }
            for j in 1..p_k+1
            {
                let mut saved = 0.0;
                for r in 0..j 
                {
                    let ri = right[r];
                    let le = left[p_k - j + r];
                    let temp = shape_ders[r] / (ri + le);
                    shape_ders[r] = saved + (ri * temp);
                    saved = le * temp;
                }
                shape_ders[j] = saved;
            }
        }

        for k2 in 1..k+1 // loop over derivative
        {
            let p_k2 = p_k + k2;
            let p_k2_f = p_k2 as f64;
            let mut saved = 0.0;
            for r in 0..p_k2
            {
                let denom = self.knots[i + 1 + r] - self.knots[i + 1 + r - p_k2];
                let temp = ( p_k2_f  / denom ) * shape_ders[r];
                shape_ders[r] = saved - temp;
                saved = temp;
            }
            shape_ders[p_k2] = saved;
        }
    }
}



// ------------------------------------------- Tests -------------------------------------------- //

#[cfg(test)]
mod tests {

    use super::BsplineBasis;
    use super::PMAX;

    use approx::assert_relative_eq;
    use serde::Deserialize;
    use std::fs;

    const ZERO_THRESHOLD: f64=1e-14;

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

    macro_rules! eval_diff {
        ($test_name:ident, $knots:ident, $ders:ident, $order:expr) => {
            #[test]
            fn $test_name() {
                let p = $order;
                let test_data = TestData::new();
                let basis = BsplineBasis::new(test_data.$knots.values.clone());
                for (idx, u) in test_data.u.values.iter().enumerate()
                {
                    let ders_start = idx * (p + 1);
                    let ders_end = (idx+1) * (p + 1);
                    let ders1_all = test_data.$ders.values[ders_start..ders_end].to_vec();

                    let (start, _end, num_basis) = basis.non_zero_basis(*u, p);

                    for j in 0..p+1 
                    {
                        let mut ders1 = ders1_all[j].clone();
                        ders1.iter_mut().for_each(|elem| {
                            if elem.abs() < ZERO_THRESHOLD {
                                *elem = 0.0;
                            }
                        });
                        let mut ders2 = [0.0; PMAX+1];
                        basis.eval_diff(*u, p, j, &mut ders2);
                        ders2.iter_mut().for_each(|elem| {
                            if elem.abs() < ZERO_THRESHOLD {
                                *elem = 0.0;
                            }
                        });
                        for k in 0..num_basis 
                        {
                            let val1 = ders1[k + start];
                            let val2 = ders2[k];
                            assert_relative_eq!(val1, val2, max_relative = 1e-14);
                        }
                    }
                }
            }
        };
    }

    eval_diff!(eval_diff0, knots_p0, ders_p0, 0);
    eval_diff!(eval_diff1, knots_p1, ders_p1, 1);
    eval_diff!(eval_diff2, knots_p2, ders_p2, 2);
    eval_diff!(eval_diff3, knots_p3, ders_p3, 3);
    eval_diff!(eval_diff4, knots_p4, ders_p4, 4);

}