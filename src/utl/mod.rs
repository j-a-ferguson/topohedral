//! This is the UTiLities module
//!
//!
//!
//!

use static_assertions::const_assert;

// -------------------------------------------- Macros ------------------------------------------ //

// ------------------------------------------- Structs ------------------------------------------ //

#[repr(C)]
pub struct LinearIndex<const N: usize> {
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

    pub fn new(dims: &[usize; N]) -> Self {
        assert!(dims.len() == N);

        let mut lin_idx = LinearIndex {
            dims: *dims,
            idx_helper: [1; N],
            tuple_helper: *dims,
        };

        for i in 1..N {
            lin_idx.idx_helper[i] = lin_idx.idx_helper[i - 1] * dims[i - 1];
            lin_idx.tuple_helper[i] = lin_idx.tuple_helper[i] * lin_idx.tuple_helper[i - 1];
        }
        lin_idx
    }

    pub fn index(&self, indices: &[usize; N]) -> usize {
        let mut idx = 0usize;
        for i in 0..N {
            idx += (indices[i] * self.idx_helper[i]);
        }
        idx
    }
    //..............................................................................................

    pub fn tuple(&self, idx: usize) -> [usize; N] {
        let mut out = [0; N];
        out[0] = idx % self.dims[0];
        for i in 1..N {
            let base = (idx / self.tuple_helper[i]) * self.tuple_helper[i];
            out[i] = (idx - base) / self.tuple_helper[i - 1];
        }
        out
    }
    //..............................................................................................
}

// ------------------------------------------- Tests -------------------------------------------- //

mod tests {
    use crate::utl::LinearIndex;

    #[test]
    fn linear_index2() {
        let lin_idx = LinearIndex::new(&[3, 4]);

        let mut idx1 = 0;
        for j in 0..4
        {
            for i in 0..3
            {
                let tuple1 = [i, j];
                let tuple2 = lin_idx.tuple(idx1);
                let idx2 = lin_idx.index(&tuple1);
                assert_eq!(idx1, idx2);
                assert_eq!(tuple1, tuple2);
                idx1 += 1;
            }
        }
    }

    #[test]
    fn linear_index3() {

        let lin_idx = LinearIndex::new(&[3, 4, 2]);

        let mut idx1 = 0;
        for k in 0..2
        {
            for j in 0..4
            {
                for i in 0..3
                {
                    let tuple1 = [i, j, k];
                    let tuple2 = lin_idx.tuple(idx1);
                    let idx2 = lin_idx.index(&tuple1);
                    assert_eq!(idx1, idx2);
                    assert_eq!(tuple1, tuple2);
                    idx1 += 1;
                }
            }
        }

    }

}
