// Copyright 2026 Jan Niklas Siemer
//
// This file is part of scalable_rbe.
//
// scalable_rbe is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! Implementation of the matrix type in coefficient embedding.

use crate::dgs::{DGS, sample_binomial};
use rand::{Rng, distr::Distribution, rngs::ThreadRng};
use rand_distr::Uniform;

/// Holds a matrix containing entries of polynomials of degree `d` in coefficent representation
/// no larger than [`u64::MAX`] with `nr_rows`-many rows and `nr_columns`-many columns.
///
/// Each polynomial is stored in `d` consecutive fields of `matrix`, i.e. `matrix[0..d]` represents the first entry.
/// The entries are stored in a column-wise order, i.e. `matrix[0..d * nr_columns]` contains all
/// entries of the first column.
///
/// Attributes:
/// - `matrix`: Contains each coefficient of the matrix
/// - `d`: Maximum degree of all polynomials in the matrix
/// - `nr_rows`: Number of rows of the matrix
/// - `nr_columns`: Number of columns of the matrix
#[derive(Clone, Debug, PartialEq)]
pub struct Mat {
    pub matrix: Vec<u64>,
    pub d: usize,
    pub nr_rows: usize,
    pub nr_columns: usize,
}

impl Mat {
    /// Generates a uniformly at random chosen matrix.
    ///
    /// Parameters:
    /// - `d`: Maximum degree of each polynomial in the matrix
    /// - `nr_rows`: Number of rows of the new matrix
    /// - `nr_columns`: Number of columns the new matrix
    /// - `dist`: [`Uniform`] distribution and range that the matrix is sampled with
    /// - `rng`: A [`ThreadRng`] to generate randomness
    pub fn sample_uniform<T: Rng>(
        d: usize,
        nr_rows: usize,
        nr_columns: usize,
        dist: &Uniform<u64>,
        rng: &mut T,
    ) -> Self {
        let len = d * nr_rows * nr_columns;

        let mut res = Vec::with_capacity(len);
        for _ in 0..len {
            res.push(dist.sample(rng));
        }

        Self {
            matrix: res,
            d,
            nr_rows,
            nr_columns,
        }
    }

    /// Generates a new matrix with values sampled from `dist` shifted by `n_half` reduced `mod_q`.
    ///
    /// Parameters:
    /// - `d`: Maximum degree of each polynomial in the matrix
    /// - `nr_rows`: Number of rows of the new matrix
    /// - `nr_columns`: Number of columns the new matrix
    /// - `dist`: The distribution that the matrix is sampled from
    /// - `n_half`: The value to shift the distribution by
    /// - `mod_q`: The modulus to reduce every entry by
    /// - `rng`: A [`ThreadRng`] to generate randomness
    pub fn sample_small<T: Rng, U: Distribution<u64>>(
        d: usize,
        nr_rows: usize,
        nr_columns: usize,
        dist: &U,
        n_half: u64,
        mod_q: u64,
        rng: &mut T,
    ) -> Self {
        let len = d * nr_rows * nr_columns;
        let shift = mod_q - n_half;

        let mut res = Vec::with_capacity(len);
        for _ in 0..len {
            res.push((shift + dist.sample(rng)) % mod_q);
        }

        Self {
            matrix: res,
            d,
            nr_rows,
            nr_columns,
        }
    }

    /// Generates a new matrix with each entry chosen discrete Gaussian.
    ///
    /// Parameters:
    /// - `d`: Maximum degree of each polynomial in the matrix
    /// - `nr_rows`: Number of rows of the new matrix
    /// - `nr_columns`: Number of columns the new matrix
    /// - `dgs`: A discrete Gaussian sampler to sample each coefficient with
    /// - `rng`: A [`ThreadRng`] to generate randomness
    pub fn sample_disc_gauss(
        d: usize,
        nr_rows: usize,
        nr_columns: usize,
        dgs: &impl DGS,
        rng: &mut ThreadRng,
    ) -> Self {
        let size = d * nr_rows * nr_columns;
        let mut vector = Vec::with_capacity(size);
        for _ in 0..size {
            vector.push(dgs.sample(rng));
        }

        Self {
            matrix: vector,
            d,
            nr_rows,
            nr_columns,
        }
    }

    /// Generates a new matrix with each entry chosen from the `0`-centered binomial `Bin(n, 0.5)`.
    ///
    /// Parameters:
    /// - `d`: Maximum degree of each polynomial in the matrix
    /// - `nr_rows`: Number of rows of the new matrix
    /// - `nr_columns`: Number of columns the new matrix
    /// - `n`: Number of trials / Binomial parameter `n`
    /// - `mod_q`: The modulus to reduce every coefficient by
    /// - `rng`: A [`ThreadRng`] to generate randomness
    pub fn sample_binomial(
        d: usize,
        nr_rows: usize,
        nr_columns: usize,
        n: usize,
        mod_q: u64,
        rng: &mut ThreadRng,
    ) -> Self {
        let mod_q: i64 = mod_q as i64;
        let size = d * nr_rows * nr_columns;
        let mut vector = Vec::with_capacity(size);
        for _ in 0..size {
            vector.push(sample_binomial(n, mod_q, rng));
        }

        Self {
            matrix: vector,
            d,
            nr_rows,
            nr_columns,
        }
    }

    /// Scalar multiplication of `self` with `scalar`.
    ///
    /// **WARNING:** We require `coefficient * scalar < q` in the single context where it's used.
    ///
    /// Parameters:
    /// - `scalar`: The scalar value to multiply `self` with
    pub fn mul_scalar(&mut self, scalar: u64) {
        for i in 0..self.matrix.len() {
            self.matrix[i] *= scalar;
        }
    }

    /// Scalar multiplication of `self` using [`f64`] with `scalar` reduced `mod_q`.
    ///
    /// **WARNING:** [`f64`] introduces floating point arithmetic and errors / noise.
    ///
    /// Parameters:
    /// - `scalar`: The scalar value to multiply `self` with
    /// - `mod_q`: The modulus
    pub fn mul_scalar_f64(&mut self, scalar: f64, mod_q: u64) {
        let q_f = mod_q as f64;
        for i in 0..self.matrix.len() {
            let x = self.matrix[i] as f64 * scalar;
            let mut r = x.round();
            r = r - (r / q_f).floor() * q_f;
            self.matrix[i] = r as u64;
        }
    }

    /// Reduces each coefficient of each entry of the matrix `mod_q`.
    ///
    /// Parameters:
    /// - `mod_q`: The modulus
    pub fn reduce(&mut self, mod_q: u64) {
        for i in 0..self.matrix.len() {
            self.matrix[i] %= mod_q;
        }
    }

    /// Adds the matrix `rhs` to `self` reduced `mod_q` reusing the memory of `self`.
    ///
    /// Parameters:
    /// - `rhs`: The other matrix to add to `self`
    /// - `mod_q`: The modulus
    ///
    /// ## Panics ...
    /// - if the maximum dimensions `d` mismatch.
    /// - if the number of rows or number of columns of the matrix mismatch.
    pub fn add_assign(&mut self, rhs: &Self, mod_q: u64) {
        assert_eq!(self.d, rhs.d);
        assert_eq!(self.nr_rows, rhs.nr_rows);
        assert_eq!(self.nr_columns, rhs.nr_columns);

        for i in 0..self.matrix.len() {
            let (r, overflow) = self.matrix[i].overflowing_add(rhs.matrix[i]);
            self.matrix[i] = r.wrapping_sub((overflow as u64) * mod_q) % mod_q;
        }
    }

    /// Subtracts the matrix `rhs` from `self` reduced `mod_q` reusing the memory of `self`.
    ///
    /// Parameters:
    /// - `rhs`: The other matrix to subtract from `self`
    /// - `mod_q`: The modulus
    ///
    /// ## Panics ...
    /// - if the maximum dimensions `d` mismatch.
    /// - if the number of rows or number of columns of the matrix mismatch.
    pub fn sub_assign(&mut self, rhs: &Self, mod_q: u64) {
        assert_eq!(self.d, rhs.d);
        assert_eq!(self.nr_rows, rhs.nr_rows);
        assert_eq!(self.nr_columns, rhs.nr_columns);

        for i in 0..self.matrix.len() {
            let a = self.matrix[i];
            let b = rhs.matrix[i];
            self.matrix[i] = if a >= b { a - b } else { a + mod_q - b };
        }
    }

    /// Returns the polynomial stored in `self` corresponding to the entry in `row` and `column`.
    ///
    /// Parameters:
    /// - `row`: The row that the entry is in
    /// - `column`: The column that the entry is in
    pub fn get_entry(&self, row: usize, column: usize) -> &[u64] {
        let index = self.d * row + self.d * self.nr_rows * column;
        &self.matrix[index..index + self.d]
    }

    /// Negates the matrix `self` reduced `mod_q`.
    ///
    /// Parameters:
    /// - `mod_q`: The modulus
    pub fn negate(&mut self, mod_q: u64) {
        for i in 0..self.matrix.len() {
            self.matrix[i] = mod_q - self.matrix[i];
        }
    }

    /// Generates the pruned Gadget matrix `G~ = I_{nr_rows} ⊗ g~^t` with `g~^t = [base^{prune_by} ... base^{k-1}]`.
    ///
    /// Parameters:
    /// - `k`: Length of the gadget vector `g~^t`
    /// - `base`: Base of the gadget matrix
    /// - `prune_by`: The number of omitted lower order entries of `g~^t`
    /// - `nr_rows`: The number of rows that `G` should have
    /// - `d`: The maximum degree of each polynomial
    ///
    /// ## Panics...
    /// - if `k <= prune_by`.
    pub fn gen_pruned_gadget(
        k: usize,
        base: u64,
        prune_by: usize,
        nr_rows: usize,
        d: usize,
    ) -> Mat {
        assert!(
            k > prune_by,
            "Width of the gadget vector `k` needs to be larger than the number of pruned columns `prune_by`."
        );

        let mut res = vec![0; nr_rows * (k - prune_by) * nr_rows * d];
        for i in 0..nr_rows {
            for j in 0..(k - prune_by) {
                res[i * (k - prune_by) * nr_rows * d + i * d + j * nr_rows * d] =
                    base.pow((j + prune_by) as u32);
            }
        }

        Mat {
            matrix: res,
            d,
            nr_rows,
            nr_columns: nr_rows * (k - prune_by),
        }
    }

    /// Generates `G~^{-1}(self)`, i.e. the inverse of the Gadget matrix `G~ = I_{nr_rows} ⊗ g~^t`
    /// with `g~^t = [base^{prune_by} ... base^{k-1}]`. This resembles the binary decomposition if `base = 2`.
    ///
    /// Parameters:
    /// - `k`: Length of the gadget vector `g~^t`
    /// - `base`: Base of the gadget matrix
    /// - `prune_by`: The number of omitted lower order entries of `g~^t`
    /// - `mod_q`: The modulus
    ///
    /// ## Panics...
    /// - if `k <= prune_by`.
    /// - if `self` isn't a column vector.
    pub fn pruned_binary_decomp(self, k: usize, base: u64, prune_by: usize, mod_q: u64) -> Mat {
        assert_eq!(self.nr_columns, 1);
        assert!(k > prune_by, "k must be > prune_by");
        assert!(base >= 2, "Base must be 2 or greater");

        let kept = k - prune_by;
        let mut out = vec![0u64; self.nr_rows * self.d * kept];

        // This is the `d_offset` vector [b/2 ... b/2] represented as a scalar
        let b_half = base / 2;

        // This is `h = G * d_offset`
        let b_k = base.pow(k as u32);

        // h_val_int = (b/2) * (b^0 + ... + b^{k-1})
        let h_val_int =
            // b_half * ( (b^k - b^0) / (b - 1) )
            b_half.wrapping_mul( (b_k.wrapping_sub(1)) / (base - 1) );

        let h_val_mod_q = h_val_int % mod_q;

        for i in 0..self.nr_rows {
            for j in 0..self.d {
                let orig_val = self.matrix[i * self.d + j];

                // Calculate (x + h) mod q
                let val_plus_h = (orig_val + h_val_mod_q) % mod_q;

                // Apply G^{-1} (flooring decomposition)
                let mut temp_val = val_plus_h;
                for _ in 0..prune_by {
                    temp_val /= base;
                }

                // Decompose the remaining 'kept' digits
                for digit in 0..kept {
                    let out_idx = i * (kept * self.d) + digit * self.d + j;

                    let decomp_digit = temp_val % base;
                    temp_val /= base;

                    // Subtract d_offset mod q
                    out[out_idx] = (decomp_digit + mod_q - b_half) % mod_q;
                }
            }
        }

        Mat {
            matrix: out,
            d: self.d,
            nr_rows: self.nr_rows * kept,
            nr_columns: 1,
        }
    }
}

#[cfg(test)]
mod test_mat {
    use super::Mat;
    use crate::ntt::MatNTT;
    use concrete_ntt::prime64::Plan;
    use rand::rng;
    use rand_distr::Uniform;

    #[test]
    fn round_trip_pruned_gadget() {
        let nr_rows = 200;
        let q = 257;
        let base = 50;
        let prune_by = 1;
        let d = 16;
        let log_q = (q as f64).log(base as f64).ceil() as usize;
        let ntt_plan = Plan::try_new(d, q).unwrap();

        let vec_y = Mat::sample_uniform(d, nr_rows, 1, &Uniform::new(0, q).unwrap(), &mut rng());

        let inv = vec_y.clone().pruned_binary_decomp(log_q, base, prune_by, q);

        let g = Mat::gen_pruned_gadget(log_q, base, prune_by, nr_rows, d);

        let mut res = MatNTT::ntt(g, &ntt_plan).mul(&MatNTT::ntt(inv, &ntt_plan), &ntt_plan);
        res.normalise(&ntt_plan);
        let res = res.inv_ntt(&ntt_plan);

        for entry_i in 0..vec_y.nr_rows {
            let orig_entry = vec_y.get_entry(entry_i, 0);
            let res_entry = res.get_entry(entry_i, 0);

            for coeff_j in 0..vec_y.d {
                let orig_coeff = orig_entry[coeff_j] as i64;
                let res_coeff = res_entry[coeff_j] as i64;

                let diff = orig_coeff - res_coeff;
                let diff = if diff < -128 { q as i64 + diff } else { diff };

                assert!(diff <= (base as i64).pow(prune_by as u32) / 2);
                assert!(diff >= -(base as i64).pow(prune_by as u32) / 2);
            }
        }
    }

    #[test]
    fn correctness_add_0() {
        let q = 17;
        let d = 4;
        let nr_rows = 2;
        let nr_columns = 3;
        let mut mat_0 = Mat {
            matrix: vec![
                1, 11, 14, 5, 6, 2, 7, 13, 6, 14, 7, 13, 2, 12, 6, 15, 16, 3, 11, 15, 15, 16, 0, 3,
            ],
            d,
            nr_rows,
            nr_columns,
        };
        let mat_1 = Mat {
            matrix: vec![
                6, 3, 5, 11, 2, 4, 6, 11, 5, 8, 4, 1, 11, 5, 8, 8, 5, 12, 7, 11, 1, 10, 15, 16,
            ],
            d,
            nr_rows,
            nr_columns,
        };
        let cmp = Mat {
            matrix: vec![
                7, 14, 2, 16, 8, 6, 13, 7, 11, 5, 11, 14, 13, 0, 14, 6, 4, 15, 1, 9, 16, 9, 15, 2,
            ],
            d,
            nr_rows,
            nr_columns,
        };

        mat_0.add_assign(&mat_1, q);

        assert_eq!(cmp, mat_0);
    }

    #[test]
    fn correctness_add_1() {
        let q = 7;
        let d = 4;
        let nr_rows = 2;
        let nr_columns = 1;
        let mut mat_0 = Mat {
            matrix: vec![1, 2, 4, 3, 6, 2, 6, 6],
            d,
            nr_rows,
            nr_columns,
        };
        let mat_1 = Mat {
            matrix: vec![2, 6, 4, 3, 2, 1, 3, 0],
            d,
            nr_rows,
            nr_columns,
        };
        let cmp = Mat {
            matrix: vec![3, 1, 1, 6, 1, 3, 2, 6],
            d,
            nr_rows,
            nr_columns,
        };

        mat_0.add_assign(&mat_1, q);

        assert_eq!(cmp, mat_0);
    }

    #[test]
    fn correctness_sub_0() {
        let q = 17;
        let d = 4;
        let nr_rows = 2;
        let nr_columns = 3;
        let mut mat_0 = Mat {
            matrix: vec![
                12, 9, 13, 5, 13, 10, 8, 4, 6, 10, 0, 5, 1, 16, 9, 10, 2, 13, 5, 10, 6, 8, 5, 7,
            ],
            d,
            nr_rows,
            nr_columns,
        };
        let mat_1 = Mat {
            matrix: vec![
                10, 3, 3, 12, 5, 10, 4, 5, 15, 16, 2, 9, 11, 7, 7, 9, 1, 15, 0, 8, 13, 7, 14, 7,
            ],
            d,
            nr_rows,
            nr_columns,
        };
        let cmp = Mat {
            matrix: vec![
                2, 6, 10, 10, 8, 0, 4, 16, 8, 11, 15, 13, 7, 9, 2, 1, 1, 15, 5, 2, 10, 1, 8, 0,
            ],
            d,
            nr_rows,
            nr_columns,
        };

        mat_0.sub_assign(&mat_1, q);

        assert_eq!(cmp, mat_0);
    }

    #[test]
    fn correctness_sub_1() {
        let q = 7;
        let d = 4;
        let nr_rows = 2;
        let nr_columns = 1;
        let mut mat_0 = Mat {
            matrix: vec![6, 1, 4, 1, 6, 0, 1, 3],
            d,
            nr_rows,
            nr_columns,
        };
        let mat_1 = Mat {
            matrix: vec![5, 1, 2, 6, 6, 6, 2, 2],
            d,
            nr_rows,
            nr_columns,
        };
        let cmp = Mat {
            matrix: vec![1, 0, 2, 2, 0, 1, 6, 1],
            d,
            nr_rows,
            nr_columns,
        };

        mat_0.sub_assign(&mat_1, q);

        assert_eq!(cmp, mat_0);
    }
}
