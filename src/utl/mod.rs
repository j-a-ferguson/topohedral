//! This is the UTiLities module
//! 
//! 
//! 
//! 

use static_assertions::const_assert;

// -------------------------------------------- Macros ------------------------------------------ //

// ------------------------------------------- Structs ------------------------------------------ //

#[repr(C)]
pub struct LinearIndex<const N: usize> 
{
    dims: [usize; N],
    idx_helper: [usize; N], 
    tuple_helper: [usize; N],
}

// ------------------------------------------- Enums -------------------------------------------- //


// ------------------------------------------- Types -------------------------------------------- //


// ----------------------------------------- Constants ------------------------------------------ //


// --------------------------------------- Free Functions---------------------------------------- //


// --------------------------------------- Implementations -------------------------------------- //

impl<const N: usize> LinearIndex<N> {

    // const_assert!(N > 0);


    pub fn new(dims: &[usize; N]) -> Self
    {
        assert!(dims.len() == N);

        let mut lin_idx = LinearIndex{
            dims: *dims,
            idx_helper: [1; N], 
            tuple_helper: *dims,
        };

        for i in 1..N
        {
            lin_idx.idx_helper[i] = lin_idx.idx_helper[i-1] * dims[i-1];
            lin_idx.tuple_helper[i] = lin_idx.tuple_helper[i] * lin_idx.tuple_helper[i-1];
        }
        lin_idx
    }

    pub fn index(&self, indices: &[usize; N])
    {
        let idx = 0usize;
        for i in 0..N
        {
            idx += (indices[i] * )
        }
    }
}

// ------------------------------------------- Tests -------------------------------------------- //