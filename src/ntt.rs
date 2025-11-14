//! Implementation of the matrix type in NTT form.

use crate::mat::Mat;
use concrete_ntt::prime64::Plan;
use rand::{
    Rng,
    distr::{Distribution, Uniform},
};

/// Holds a matrix in NTT form containing entries of polynomials of degree `d`
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
pub struct MatNTT {
    pub matrix: Vec<u64>,
    pub d: usize,
    pub nr_rows: usize,
    pub nr_columns: usize,
}

impl MatNTT {
    /// Computes the transform from coefficient representation to NTT form.
    ///
    /// **WARNING:** On entry, the buffer holds the polynomial coefficients in standard order.
    /// On exit, the buffer holds the negacyclic NTT transform coefficients in bit reversed order.
    ///
    /// Parameters:
    /// - `matrix`: The matrix to transform to NTT form
    /// - `ntt_plan`: The NTT-[`Plan`] that stores several information about NTT transform, etc.
    ///
    /// ## Panics...
    /// - if the size of the polynomials mismatches the polynomials assumed in `ntt_plan`.
    /// - if a value larger than `q` of `ntt_plan` occurs.
    pub fn ntt(mut matrix: Mat, ntt_plan: &Plan) -> Self {
        for chunk in (0..matrix.matrix.len()).step_by(matrix.d) {
            ntt_plan.fwd(&mut matrix.matrix[chunk..chunk + matrix.d]);
        }
        // Only required if no action is performed on it, which should never be the case
        // ntt_plan.normalize(&mut res);

        MatNTT {
            matrix: matrix.matrix,
            d: matrix.d,
            nr_rows: matrix.nr_rows,
            nr_columns: matrix.nr_columns,
        }
    }

    /// Computes the inverse NTT transform, i.e. moves `self` from NTT form to coefficient representation.
    ///
    /// **WARNING:** On entry, the buffer holds the negacyclic NTT transform coefficients in bit reversed order.
    /// On exit, the buffer holds the polynomial coefficients in standard order.
    /// Also, it's required that the input is `normalised` / reduced modulo the appropriate `q`.
    ///
    /// Parameters:
    /// - `ntt_plan`: The NTT-[`Plan`] that stores several information about NTT transform, etc.
    ///
    /// ## Panics...
    /// - if the size of the polynomials mismatches the polynomials assumed in `ntt_plan`.
    /// - if a value larger than `q` of `ntt_plan` occurs.
    pub fn inv_ntt(mut self, ntt_plan: &Plan) -> Mat {
        for chunk in (0..self.matrix.len()).step_by(self.d) {
            ntt_plan.inv(&mut self.matrix[chunk..chunk + self.d]);
        }

        Mat {
            matrix: self.matrix,
            d: self.d,
            nr_rows: self.nr_rows,
            nr_columns: self.nr_columns,
        }
    }

    /// Multiplies the matrix `self` to `rhs` using NTT / component-wise multiplication.
    ///
    /// Parameters:
    /// - `rhs`: The other matrix to multiply from the right-hand side with
    /// - `ntt_plan`: The NTT-[`Plan`] that stores several information about NTT transform, etc.
    ///
    /// ## Panics...
    /// - if the size of the polynomials mismatches.
    /// - if the number of columns of `self` mismatches the number of rows of `rhs`.
    pub fn mul(&self, rhs: &MatNTT, ntt_plan: &Plan) -> MatNTT {
        assert_eq!(self.d, rhs.d);
        assert_eq!(self.nr_columns, rhs.nr_rows);

        let mut res = MatNTT {
            matrix: vec![0; self.d * self.nr_rows * rhs.nr_columns],
            d: self.d,
            nr_rows: self.nr_rows,
            nr_columns: rhs.nr_columns,
        };

        for column in 0..rhs.nr_columns {
            for row in 0..self.nr_rows {
                for i in 0..self.nr_columns {
                    ntt_plan.mul_accumulate(
                        res.get_entry_mut(row, column),
                        self.get_entry(row, i),
                        rhs.get_entry(i, column),
                    );
                }
            }
        }

        res
    }

    /// Adds the matrix `self` to `rhs` reusing the memory of `self`.
    ///
    /// **WARNING:** This function assumes that for every entry `self.matrix[i] + rhs.matrix[i] < u64::MAX`.
    /// Furthermore, it does not reduce mod `q`.
    ///
    /// Parameters:
    /// - `rhs`: The other matrix to add to `self`
    ///
    /// ## Panics...
    /// - if the size of the polynomials mismatches.
    /// - if the number of rows or the number of columns of `self` and `rhs` mismatch.
    pub fn add_assign(&mut self, rhs: &MatNTT) {
        assert_eq!(self.d, rhs.d);
        assert_eq!(self.nr_rows, rhs.nr_rows);
        assert_eq!(self.nr_columns, rhs.nr_columns);

        for i in 0..self.matrix.len() {
            self.matrix[i] += rhs.matrix[i];
        }
    }

    /// Subtracts the matrix `rhs` from `self` reusing the memory of `self`.
    ///
    /// Parameters:
    /// - `rhs`: The other matrix to subtract from `self`
    /// - `mod_q`: The modulus to reduce by
    ///
    /// ## Panics...
    /// - if the size of the polynomials mismatches.
    /// - if the number of rows or the number of columns of `self` and `rhs` mismatch.
    pub fn sub_assign(&mut self, rhs: &MatNTT, mod_q: u64) {
        assert_eq!(self.d, rhs.d);
        assert_eq!(self.nr_rows, rhs.nr_rows);
        assert_eq!(self.nr_columns, rhs.nr_columns);

        for i in 0..self.matrix.len() {
            let a = self.matrix[i];
            let b = rhs.matrix[i];
            self.matrix[i] = if a >= b { a - b } else { a + mod_q - b };
        }
    }

    /// Generates a matrix in NTT form with each entry sampled uniformly at random.
    ///
    /// Parameters:
    /// - `d`: Maximum degree of each polynomial in the matrix
    /// - `nr_rows`: Number of rows of the new matrix
    /// - `nr_columns`: Number of columns the new matrix
    /// - `mod_q`: The modulus that defines the size of the interval to sample from
    /// - `rng`: A [`Rng`] to generate randomness
    pub fn sample_uniform<T: Rng>(
        d: usize,
        nr_rows: usize,
        nr_columns: usize,
        mod_q: u64,
        rng: &mut T,
    ) -> MatNTT {
        let len = d * nr_rows * nr_columns;
        let dist = Uniform::new(0, mod_q).unwrap();

        let mut res = Vec::with_capacity(len);
        for _ in 0..len {
            res.push(dist.sample(rng));
        }

        MatNTT {
            matrix: res,
            d,
            nr_rows,
            nr_columns,
        }
    }

    /// Returns the polynomial stored in `self` corresponding to the entry in `row` and `column`
    /// as a mutable pointer.
    ///
    /// Parameters:
    /// - `row`: The row that the entry is in
    /// - `column`: The column that the entry is in
    fn get_entry_mut(&mut self, row: usize, column: usize) -> &mut [u64] {
        let index = self.d * row + self.d * self.nr_rows * column;
        &mut self.matrix[index..index + self.d]
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

    /// Normalises each polynomial, i.e. reduces it mod `q`.
    ///
    /// Parameters:
    /// - `ntt_plan`: The NTT-[`Plan`] that stores several information about NTT transform, etc.
    pub fn normalise(&mut self, ntt_plan: &Plan) {
        ntt_plan.normalize(&mut self.matrix);
    }
}

#[cfg(test)]
mod test {
    use crate::ntt::MatNTT;
    use concrete_ntt::prime64::Plan;

    #[test]
    fn correctness_mul_0() {
        let q = 257;
        let d = 16;
        let ntt_plan = Plan::try_new(d as usize, q).unwrap();
        let mat_ntt_0 = MatNTT {
            matrix: vec![
                84, 214, 240, 52, 27, 68, 185, 256, 104, 242, 131, 183, 154, 78, 34, 203, 146, 113,
                5, 60, 72, 3, 168, 229, 227, 223, 94, 192, 115, 233, 86, 192,
            ],
            d,
            nr_rows: 1,
            nr_columns: 2,
        };
        let mat_ntt_1 = MatNTT {
            matrix: vec![
                105, 89, 245, 111, 140, 142, 221, 17, 211, 125, 243, 27, 44, 214, 41, 187, 166, 78,
                239, 36, 251, 48, 194, 14, 159, 62, 175, 19, 185, 251, 44, 92,
            ],
            d,
            nr_rows: 2,
            nr_columns: 1,
        };
        let cmp = MatNTT {
            matrix: vec![
                10, 135, 232, 46, 145, 227, 143, 119, 206, 249, 14, 71, 163, 217, 163, 248,
            ],
            d,
            nr_rows: 1,
            nr_columns: 1,
        };

        let mut res_ntt = mat_ntt_0.mul(&mat_ntt_1, &ntt_plan);
        res_ntt.normalise(&ntt_plan);

        assert_eq!(cmp, res_ntt);
    }

    #[test]
    fn correctness_mul_1() {
        let q = 257;
        let d = 16;
        let ntt_plan = Plan::try_new(d as usize, q).unwrap();
        let mat_ntt_0 = MatNTT {
            matrix: vec![
                176, 87, 175, 206, 3, 8, 2, 90, 33, 123, 247, 116, 34, 15, 134, 203, 106, 56, 5,
                97, 183, 16, 34, 132, 0, 110, 202, 143, 128, 13, 119, 219,
            ],
            d,
            nr_rows: 2,
            nr_columns: 1,
        };
        let mat_ntt_1 = MatNTT {
            matrix: vec![
                71, 220, 240, 8, 154, 146, 184, 198, 207, 210, 180, 208, 135, 203, 139, 89, 95, 82,
                191, 23, 246, 130, 195, 210, 219, 248, 5, 73, 218, 234, 44, 1,
            ],
            d,
            nr_rows: 1,
            nr_columns: 2,
        };
        let cmp = MatNTT {
            matrix: vec![
                10, 104, 55, 103, 61, 73, 23, 150, 186, 233, 16, 223, 62, 110, 104, 53, 117, 256,
                75, 177, 123, 146, 134, 220, 0, 223, 88, 60, 52, 181, 54, 142, 17, 221, 17, 7, 14,
                65, 185, 89, 18, 236, 29, 208, 142, 123, 240, 93, 19, 30, 140, 27, 83, 130, 61, 62,
                0, 163, 31, 26, 202, 158, 6, 94,
            ],
            d,
            nr_rows: 2,
            nr_columns: 2,
        };

        let mut res_ntt = mat_ntt_0.mul(&mat_ntt_1, &ntt_plan);
        res_ntt.normalise(&ntt_plan);

        assert_eq!(cmp, res_ntt);
    }

    #[test]
    fn correctness_add_0() {
        let q = 257;
        let d = 16;
        let ntt_plan = Plan::try_new(d, q).unwrap();
        let mut mat_ntt_0 = MatNTT {
            matrix: vec![
                165, 194, 102, 136, 178, 97, 236, 119, 91, 235, 150, 169, 138, 25, 226, 80, 200,
                108, 53, 5, 101, 244, 29, 107, 62, 160, 239, 220, 93, 174, 26, 5,
            ],
            d,
            nr_rows: 2,
            nr_columns: 1,
        };
        let mat_ntt_1 = MatNTT {
            matrix: vec![
                160, 108, 29, 64, 21, 217, 221, 208, 57, 122, 250, 154, 158, 256, 197, 143, 148,
                192, 99, 199, 140, 115, 72, 25, 82, 8, 77, 44, 98, 190, 36, 254,
            ],
            d,
            nr_rows: 2,
            nr_columns: 1,
        };
        let cmp = MatNTT {
            matrix: vec![
                197, 51, 217, 141, 157, 116, 141, 165, 202, 199, 25, 229, 147, 130, 171, 30, 86,
                83, 138, 77, 256, 167, 183, 201, 9, 139, 84, 145, 28, 87, 36, 225,
            ],
            d,
            nr_rows: 2,
            nr_columns: 1,
        };

        mat_ntt_0.add_assign(&mat_ntt_1);
        mat_ntt_0.normalise(&ntt_plan);

        assert_eq!(cmp, mat_ntt_0);
    }

    #[test]
    fn correctness_add_1() {
        let q = 257;
        let d = 16;
        let ntt_plan = Plan::try_new(d, q).unwrap();
        let mut mat_ntt_0 = MatNTT {
            matrix: vec![
                150, 95, 252, 171, 221, 149, 198, 12, 201, 218, 62, 206, 14, 178, 24, 256, 142, 7,
                61, 205, 196, 232, 134, 237, 170, 37, 197, 95, 223, 13, 106, 25,
            ],
            d,
            nr_rows: 1,
            nr_columns: 2,
        };
        let mat_ntt_1 = MatNTT {
            matrix: vec![
                74, 37, 144, 241, 251, 184, 206, 172, 24, 91, 30, 70, 66, 234, 158, 38, 36, 94, 14,
                221, 144, 91, 34, 159, 207, 248, 93, 119, 142, 155, 173, 253,
            ],
            d,
            nr_rows: 1,
            nr_columns: 2,
        };
        let cmp = MatNTT {
            matrix: vec![
                14, 201, 89, 90, 158, 69, 218, 140, 255, 196, 70, 210, 5, 90, 172, 179, 236, 183,
                85, 123, 214, 229, 139, 89, 136, 66, 243, 174, 71, 139, 162, 178,
            ],
            d,
            nr_rows: 1,
            nr_columns: 2,
        };

        mat_ntt_0.add_assign(&mat_ntt_1);
        mat_ntt_0.normalise(&ntt_plan);

        assert_eq!(cmp, mat_ntt_0);
    }

    #[test]
    fn correctness_sub_0() {
        let q = 257;
        let d = 16;
        let ntt_plan = Plan::try_new(d, q).unwrap();
        let mut mat_ntt_0 = MatNTT {
            matrix: vec![
                9, 110, 36, 226, 252, 67, 22, 37, 235, 124, 12, 10, 28, 219, 183, 23, 117, 132,
                223, 31, 8, 255, 223, 216, 197, 253, 89, 30, 89, 201, 144, 79,
            ],
            d,
            nr_rows: 2,
            nr_columns: 1,
        };
        let mat_ntt_1 = MatNTT {
            matrix: vec![
                133, 63, 77, 82, 201, 210, 42, 89, 137, 29, 137, 249, 41, 99, 29, 112, 213, 61,
                205, 19, 209, 187, 56, 85, 35, 122, 58, 168, 113, 188, 30, 6,
            ],
            d,
            nr_rows: 2,
            nr_columns: 1,
        };
        let cmp = MatNTT {
            matrix: vec![
                185, 19, 142, 9, 212, 232, 63, 61, 231, 22, 201, 226, 208, 136, 106, 139, 251, 149,
                226, 65, 132, 197, 155, 217, 235, 217, 18, 152, 127, 49, 232, 117,
            ],
            d,
            nr_rows: 2,
            nr_columns: 1,
        };

        mat_ntt_0.sub_assign(&mat_ntt_1, q);
        mat_ntt_0.normalise(&ntt_plan);

        assert_eq!(cmp, mat_ntt_0);
    }

    #[test]
    fn correctness_sub_1() {
        let q = 257;
        let d = 16;
        let ntt_plan = Plan::try_new(d, q).unwrap();
        let mut mat_ntt_0 = MatNTT {
            matrix: vec![
                64, 64, 145, 135, 98, 156, 86, 250, 179, 239, 201, 193, 184, 231, 2, 185, 83, 135,
                193, 144, 85, 62, 128, 120, 18, 204, 194, 68, 58, 190, 32, 89,
            ],
            d,
            nr_rows: 1,
            nr_columns: 2,
        };
        let mat_ntt_1 = MatNTT {
            matrix: vec![
                98, 85, 233, 57, 118, 140, 144, 120, 116, 44, 132, 71, 17, 215, 12, 14, 95, 161,
                96, 192, 13, 140, 205, 163, 20, 208, 239, 239, 147, 219, 57, 230,
            ],
            d,
            nr_rows: 1,
            nr_columns: 2,
        };
        let cmp = MatNTT {
            matrix: vec![
                30, 79, 123, 37, 63, 1, 157, 233, 20, 221, 181, 104, 155, 1, 160, 91, 192, 159,
                247, 254, 133, 220, 204, 174, 32, 64, 206, 166, 139, 207, 143, 200,
            ],
            d,
            nr_rows: 1,
            nr_columns: 2,
        };

        mat_ntt_0.sub_assign(&mat_ntt_1, q);
        mat_ntt_0.normalise(&ntt_plan);

        assert_eq!(cmp, mat_ntt_0);
    }
}
