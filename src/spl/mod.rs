
use float_cmp::approx_eq;

const PMAX: usize = 8;

pub struct BsplineBasis
{
    knots: Vec<f64>,
}


impl BsplineBasis
{
    pub fn new<I: IntoIterator<Item=f64>>(start: I) -> Self 
    {
        let knots: Vec<f64> = start.into_iter().collect();
        Self{knots}
    }

    pub fn eval(&self, u: f64, i: usize, p: usize) -> f64
    {
        let umin = self.knots.first().unwrap();
        let umax = self.knots.last().unwrap();
        assert![!self.knots.is_empty()];
        assert![u >= *umin && u <= *umax];
        assert![p < PMAX];

        let mut basis_fun_i = 0.0;
        let n = self.knots.len();

        let u_is_u0 = approx_eq!(f64, u, *umin, ulps = 4); 
        let u_is_un =  approx_eq!(f64, u, *umax, ulps = 4); 

        if (u_is_u0 && i == 0) || (u_is_un && i == (n - p - 2)){
            basis_fun_i = 1.0;
        }
        else {

            let mut njp: [f64; PMAX] = [0.0; PMAX];

            for j in (0usize..((p+1)as usize))
            {
                njp[j] = 1.0
            }

            let (mut saved, mut temp, mut uleft, mut uright) = (0.0, 0.0, 0.0, 0.0);

            for k in 1..p+1
            {
                if njp[0] == 0.0
                {
                    saved = 0.0
                }
                else 
                {
                    saved = ((u - self.knots[i]) * njp[0]) / (self.knots[i + k] - self.knots[i])
                }

                for j in 0..(p - k + 1)
                {
                    uleft = self.knots[i + j + 1];
                    uright = self.knots[i + j + k + 1];
                    if njp[j + 1] == 0.0
                    {
                        njp[j] = saved;
                        saved = 0.0;
                    }
                    else {
                        temp = njp[j + 1] / (uright - uleft);
                        njp[j] = saved + (uright - u) * temp;
                        saved = (u - uleft) * temp;
                    }
                }
            }
            basis_fun_i = njp[0];
        }
        return basis_fun_i;
    }
}



#[cfg(test)]
mod tests {
    use super::BsplineBasis;

    #[test]
    fn construction()
    {
        let knots = vec![0.0, 0.1, 0.2, 0.3];
        let bspline_basis = BsplineBasis::new(knots.clone());
        assert_eq!(bspline_basis.knots, knots);

    }

    #[test]
    fn eval
}