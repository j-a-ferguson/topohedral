//! This is the UTiLities module
//!
//!
//!
//!

use static_assertions::const_assert;
use std::ops::{Index, IndexMut};
use std::fmt::{self, write};

// -------------------------------------------- Macros ------------------------------------------ //

// ------------------------------------------- Structs ------------------------------------------ //

#[repr(C)]
pub struct NDArray<'a, T, const N: usize> {
    data: &'a mut [T],
    dims: [usize; N],
    idx_helper: [usize; N],
    tuple_helper: [usize; N],
}

// ------------------------------------------- Enums -------------------------------------------- //

// ------------------------------------------- Types -------------------------------------------- //

// ----------------------------------------- Constants ------------------------------------------ //

const MAX_DIM: usize = 4;

// --------------------------------------- Free Functions---------------------------------------- //

// --------------------------------------- Implementations -------------------------------------- //

impl<'a, T, const N: usize> NDArray<'a, T, N> {
    // const_assert!(N > 0);

    pub fn new(data: &'a mut [T], dims: &[usize]) -> Self {
        debug_assert!(dims.len() >= N);
        debug_assert!(data.len() >= dims.iter().product());

        let mut nd_array = NDArray::<'a, T, N>{
            data: data,
            dims: [0; N],
            idx_helper: [1; N],
            tuple_helper: [0; N],
        };

        nd_array.dims.clone_from_slice(&dims[0..N]);
        nd_array.tuple_helper.clone_from_slice(&dims[0..N]);

        for i in 1..N {
            nd_array.idx_helper[i] = nd_array.idx_helper[i - 1] * dims[i - 1];
            nd_array.tuple_helper[i] = nd_array.tuple_helper[i] * nd_array.tuple_helper[i - 1];
        }
        nd_array
    }

    pub fn lin_index(&self, indices: &[usize]) -> usize {
        debug_assert!(indices.len() >= N);
        let mut idx = 0usize;
        for i in 0..N {
            idx += (indices[i] * self.idx_helper[i]);
        }
        idx
    }
    //..............................................................................................

    pub fn tuple_index(&self, idx: usize) -> [usize; N] {
        let mut out = [0; N];
        out[0] = idx % self.dims[0];
        for i in 1..N {
            let base = (idx / self.tuple_helper[i]) * self.tuple_helper[i];
            out[i] = (idx - base) / self.tuple_helper[i - 1];
        }
        out
    }
    //.............................................................................................
}

impl<'a, T, const N: usize> Index<&[usize; N]> for NDArray<'a, T, N>
{
    type Output = T;

    fn index(&self, index_tuple: &[usize; N]) -> &Self::Output {
        let idx = self.lin_index(index_tuple);
        &self.data[idx]
    }
}

impl<'a, T, const N: usize> IndexMut<&[usize; N]> for NDArray<'a, T, N>
{
    fn index_mut(&mut self, index_tuple: &[usize; N]) -> &mut Self::Output {
        let idx = self.lin_index(index_tuple);
        &mut self.data[idx]
    }
}

impl<'a, T, const N: usize> fmt::Display for NDArray<'a, T, N>
where 
    T: fmt::Display
{

    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result
    {
        let len = self.dims.iter().product::<usize>();

        match N {
            1 => {
                write!(f, "[")?;
                for i in 0..len-1
                {
                    let val = &self.data[i];
                    write!(f, "{}, ", val)?;
                }
                let val = self.data.last().unwrap();
                write!(f, "{}]", val)?;
            },
            2 => {

                writeln!(f, "[")?;
                for i in 0..self.dims[0]
                {
                    write!(f, " [")?;
                    for j in 0.. self.dims[1]-1
                    {
                        let idx = self.lin_index(&[i, j]);
                        let val = &self.data[idx];
                        write!(f, "{}, ", val)?;
                    }
                    let idx = self.lin_index(&[i, self.dims[1]-1]);
                    let val = &self.data[idx];
                    writeln!(f, "{}]", val)?;
                }
                writeln!(f, "]")?;
            }
            _ if N > 2 => {

                for i in 0..len
                {
                    let i2 = self.tuple_index(i);
                    writeln!(f, "{:?}: {}", i2, self.index(&i2))?;
                }
            },
            _ => {
                panic!("N should not be 0!")
            }
        }

        Ok(())
    }
}
// ------------------------------------------- Tests -------------------------------------------- //

mod tests {
    use crate::utl::NDArray;

    #[test]
    fn linear_index2() {

        let mut data: Vec<f64> = (0..12).map(|n| n as f64).collect();
        let lin_idx = NDArray::new(data.as_mut_slice(), &[3, 4]);

        let mut idx1 = 0;
        let mut val1 = 0.0;
        for j in 0..4
        {
            for i in 0..3
            {
                let tuple1 = [i, j];
                let tuple2 = lin_idx.tuple_index(idx1);
                let idx2 = lin_idx.lin_index(&tuple1);
                let val2 = lin_idx[&tuple1];
                assert_eq!(idx1, idx2);
                assert_eq!(tuple1, tuple2);
                assert_eq!(val1, val2);
                idx1 += 1;
                val1 += 1.0;
            }
        }
    }

    #[test]
    fn linear_index3() {

        let mut data: Vec<f64> = (0..24).map(|n| n as f64).collect();
        let lin_idx = NDArray::new(data.as_mut_slice(), &[3, 4, 2]);

        let mut idx1 = 0;
        let mut val1 = 0.0;
        for k in 0..2
        {
            for j in 0..4
            {
                for i in 0..3
                {
                    let tuple1 = [i, j, k];
                    let tuple2 = lin_idx.tuple_index(idx1);
                    let idx2 = lin_idx.lin_index(&tuple1);
                    let val2 = lin_idx[&tuple1];
                    assert_eq!(idx1, idx2);
                    assert_eq!(tuple1, tuple2);
                    assert_eq!(val1, val2);
                    idx1 += 1;
                    val1 += 1.0;
                }
            }
        }

    } 



    #[test]
    fn print1() {
        let mut data: Vec<f64> = (0..24).map(|n| n as f64).collect();
        let lin_idx: NDArray<'_, f64, 1> = NDArray::new(data.as_mut_slice(), &[24]);
        println!("{}", lin_idx);
    }

    #[test]
    fn print2() {
        let mut data: Vec<f64> = (0..24).map(|n| n as f64).collect();
        let lin_idx: NDArray<'_, f64, 2> = NDArray::new(data.as_mut_slice(), &[4, 6]);
        println!("{}", lin_idx);
    }

    #[test]
    fn print3() {
        let mut data: Vec<f64> = (0..24).map(|n| n as f64).collect();
        let lin_idx: NDArray<'_, f64, 3> = NDArray::new(data.as_mut_slice(), &[3, 4, 2]);
        println!("{}", lin_idx);
    }

}
