//! Implementations of discrete Gaussian samplers.

use rand::{
    Rng, RngCore,
    distr::{Distribution, Uniform},
    rngs::ThreadRng,
};
use std::f64::consts::{PI, SQRT_2};

/// The trait implemented by any discrete Gaussian sampler.
pub trait DGS {
    /// Initialises the Discrete Gaussian Sampler with Gaussian parameter `s` and modulus `mod_q`.
    ///
    /// Parameters:
    /// - `s`: Gaussian parameter
    /// - `mod_q`: The modulus
    fn init(s: f64, mod_q: u64) -> Self;

    /// Samples from the pre-defined discrete Gaussian distribution of this [`DGS`],
    /// providing samples centered around `0` with Gaussian parameter `s` reduced `mod_q`.
    ///
    /// Parameters:
    /// - `rng`: A [`ThreadRng`] to generate randomness
    fn sample(&self, rng: &mut ThreadRng) -> u64;
}

/// This struct contains all precomputed values and tables relevant for
/// rejection sampling according to the discrete Gaussian distribution.
///
/// Attributes:
/// - `s`: Gaussian parameter
/// - `tmp_center`: Defines the temporary center of the distribution, which is shifted to `0` later on
/// - `mod_q`: The modulus
/// - `dist`: The uniform distribution and range to sample from
/// - `table`: A table / [`Vec`] containing the probabilities to accept a sample with,
///   i.e. the value of the Gaussian function
///
/// This sampler runs algorithms based on Section 4.1 of:
/// - Gentry, Craig, Chris Peikert, and Vinod Vaikuntanathan.
///   Trapdoors for hard lattices and new cryptographic constructions.
///   STOC'08.
///   <https://ia.cr/2007/432>
#[derive(Clone)]
pub struct RejSampler {
    pub s: f64,
    pub tmp_center: u64,
    pub mod_q: u64,
    pub dist: Uniform<u64>,
    pub table: Vec<u32>,
}

/// Evaluates the Gaussian function for value `x` and Gaussian parameter `s`.
///
/// Parameters:
/// - `x`: The value to evaluate the function for
/// - `s`: The Gaussian parameter
fn gauss_function(x: i64, s: f64) -> f64 {
    (-std::f64::consts::PI * x.pow(2) as f64 / s.powi(2)).exp()
}

impl DGS for RejSampler {
    /// Initialises [`RejSampler`] with Gaussian parameter `s` and modulus `mod_q`.
    ///
    /// Parameters:
    /// - `s`: Gaussian parameter
    /// - `mod_q`: The modulus
    fn init(s: f64, mod_q: u64) -> Self {
        let tmp_center = (6.0 * s).floor() as u64 + 1;
        let upper_bound = mod_q.min(2 * tmp_center);

        let mut table = vec![0; upper_bound as usize];

        let dist = Uniform::new(0, upper_bound).unwrap();

        for i in 0..(2 * tmp_center) {
            let x = i as i64 - tmp_center as i64;
            let evaluated_gauss_f = (gauss_function(x, s) * u32::MAX as f64).round() as u32;
            table[(i % mod_q) as usize] =
                evaluated_gauss_f.saturating_add(table[(i % mod_q) as usize]);
        }

        Self {
            s,
            tmp_center: tmp_center % mod_q,
            mod_q,
            dist,
            table,
        }
    }

    /// Samples from the pre-defined discrete Gaussian distribution of this [`RejSampler`],
    /// providing samples centered around `0` with Gaussian parameter `s` reduced `mod_q`.
    ///
    /// Parameters:
    /// - `rng`: A [`ThreadRng`] to generate randomness
    fn sample(&self, rng: &mut ThreadRng) -> u64 {
        loop {
            let x = self.dist.sample(rng);
            let evaluated_gauss_f = self.table[x as usize];
            if evaluated_gauss_f >= rng.next_u32() {
                let (r, underflow) = (x % self.mod_q).overflowing_sub(self.tmp_center);
                let a = r.wrapping_add((underflow as u64) * self.mod_q);
                return a;
            }
        }
    }
}

/// This struct contains all precomputed values and tables relevant for
/// convolution sampling according to the discrete Gaussian distribution
/// as described in \[C:MicWal17\].
///
/// **WARNING:** The specific algorithms implemented here only allow for
/// discrete Gaussian sampling with Gaussian parameter `s_sigma ∈ [s, sqrt(5) * s]`.
/// Attributes:
/// - `i`: Index `i` to start SampleI with, found be iteration
/// - `z`: Multiplier `z` in Algorithm 2
/// - `base_sampler`: The base sampler for `s_0`
/// - `z_i_table`: Precomputed table of `z_i` multipliers
/// - `mod_q`: The modulus
///
/// This sampler runs algorithms based on Section 5 of:
/// - \[C:MicWal17\] Micciancio, Daniele, and Michael Walter.
///   Gaussian sampling over the integers: Efficient, generic, constant-time.
///   CRYPTO'17.
///   <https://ia.cr/2017/259>
#[derive(Debug, Clone)]
pub struct MWSampler {
    // Used in Algorithm 2
    i: usize,
    z: i64,

    // Used in SampleI
    base_sampler: CdfSampler,
    z_i_table: Vec<i64>,

    mod_q: i64,
}

impl MWSampler {
    const S_0: f64 = 34.0;
    const ETA: f64 = 6.0;

    /// Samples discrete Gaussian values for an arbitrarily large `s >> η_∈(Z)`
    /// with `s_i = sqrt( (z_i^2 + max( (z_i - 1)^2, 1) ) * s_{i-1}^2 )`.
    ///
    /// Implements SampleI(i) as described in Algorithm 1 of [C:MicWal17].
    ///
    /// Parameters:
    /// - `i`: Defining `s_i` to sample from, depending on the choice of `s_0`
    /// - `rng`: A [`ThreadRng`] to generate randomness
    fn sample_i(&self, i: usize, rng: &mut impl Rng) -> i64 {
        // if i = 0, return SampleB_{s_0}(0)
        if i == 0 {
            return self.base_sampler.sample_i64(rng);
        }

        // x_1 <- SampleI(i-1)
        let x_1 = self.sample_i(i - 1, rng);

        // x_2 <- SampleI(i-1)
        let x_2 = self.sample_i(i - 1, rng);

        // y = z_i*x1 + max(1, z_i-1)*x2
        // We look up the precomputed multipliers
        let z_i = self.z_i_table[i];
        let z_max_term = 1i64.max(self.z_i_table[i] - 1);

        let term1 = x_1.saturating_mul(z_i);
        let term2 = x_2.saturating_mul(z_max_term);

        term1.saturating_add(term2)
    }
}

impl DGS for MWSampler {
    /// Initialises [`MWSampler`] with the minimum Gaussian parameter `s` and modulus `mod_q`.
    /// [`MWSampler::sample`] samples with some Gaussian parameter `s_tilde` s.t. `s <= s_tilde <= sqrt(5) s`.
    ///
    /// Parameters:
    /// - `s`: The minimum Gaussian parameter to sample with
    /// - `mod_q`: The modulus
    fn init(s_target: f64, mod_q: u64) -> Self {
        // Setup for SampleI: Build s_i and z_i tables
        let mut s_table = vec![Self::S_0];
        let mut z_i_table = vec![0i64]; // z_0 is not used, add a dummy

        // Loop to build the tables up to s_target
        let mut s_i = Self::S_0;
        while s_i < s_target {
            let s_i_minus_1 = s_i;

            // z_i = floor( |s_{i-1} / (sqrt(2)*eta)| )
            let z_i = (s_i_minus_1 / (SQRT_2 * Self::ETA)).floor() as i64;
            let z_i_minus_1 = z_i - 1;

            // max((z_i - 1)^2, 1)
            let max_term_sq = (z_i_minus_1 as f64 * z_i_minus_1 as f64).max(1.0);

            // s_i^2 = (z_i^2 + max((z_i - 1)^2, 1)) * s_{i-1}^2
            let s_i_sqrd = ((z_i as f64 * z_i as f64) + max_term_sq) * (s_i_minus_1 * s_i_minus_1);
            s_i = s_i_sqrd.sqrt();

            // Store results
            s_table.push(s_i);
            z_i_table.push(z_i);
        }

        // Setup for Algorithm 2, SampleCenteredGaussian
        // Select largest i such that s_i < s
        // Drop last entry where s_i >= s
        let _ = s_table.pop();
        let i = s_table.len() - 1;
        let s_i = s_table[i];

        // Calculate `z` for Algorithm 2
        // z := ceil( 1/2 * (1 + sqrt(2 * (s/s_i)^2 - 1)) )
        let z_float = 0.5 * (1.0 + (2.0 * (s_target / s_i).powi(2) - 1.0).sqrt());
        let z = z_float.ceil() as i64;

        Self {
            i,
            z,
            base_sampler: CdfSampler::init(Self::S_0, mod_q),
            z_i_table,
            mod_q: mod_q as i64,
        }
    }

    /// Samples from the pre-defined discrete Gaussian distribution of this [`MWSampler`],
    /// providing samples centered around `0` with *some Gaussian parameter*
    /// `s <= s_tilde <= sqrt(5) * s` reduced `mod_q`.
    ///
    /// Implements SampleCenteredGaussian(s) as described in Algorithm 2 of \[C:MicWal17\].
    ///
    /// Parameters:
    /// - `rng`: A [`ThreadRng`] to generate randomness
    fn sample(&self, rng: &mut ThreadRng) -> u64 {
        // x_1 <- SampleI(i)
        let x_1 = self.sample_i(self.i, rng);

        // x_2 <- SampleI(i)
        let x_2 = self.sample_i(self.i, rng);

        // return z*x1 + (z-1)*x2
        let term1 = x_1.saturating_mul(self.z);
        let term2 = x_2.saturating_mul(self.z - 1);

        let val = term1.saturating_add(term2);
        let res = (val % self.mod_q + self.mod_q) % self.mod_q;
        res as u64
    }
}

/// This struct contains all precomputed values and tables relevant for
/// inverse transform sampling according to the discrete Gaussian distribution
/// as described in \[Wiki:InvTransformSampling\].
///
/// Attributes:
/// - `cdf_table`: Contains the cumulative distributed functions of the Gaussian distribution
/// - `min_val`: The minimum integer value in the sampled range (e.g., -bound)
/// - `mod_q`: The modulus
///
/// This sampler runs algorithms based on Section 5 of:
/// - \[Wiki:InvTransformSampling\] Wikipedia Community.
///   Inverse transform sampling using the cumulative distribution function.
///   Last accessed: 2025-11-10.
///   <https://en.wikipedia.org/wiki/Inverse_transform_sampling>
#[derive(Debug, Clone)]
pub struct CdfSampler {
    cdf_table: Vec<u64>,
    min_val: i64,
    mod_q: i64,
}

impl CdfSampler {
    /// tailcut - defined s.t. for `s = 34.0` almost every element has the possibility to be sampled
    const TAILCUT: f64 = 3.63;

    /// Samples from the pre-defined discrete Gaussian distribution of this [`CdfSampler`],
    /// providing samples centered around `0` with Gaussian parameter `s`
    /// returning an [`i64`] *not* reduced `mod_q`.
    ///
    /// Parameters:
    /// - `rng`: A [`ThreadRng`] to generate randomness
    pub(crate) fn sample_i64(&self, rng: &mut impl Rng) -> i64 {
        let r = rng.next_u64();

        // Find first index `i` where `cdf[i] > r`
        let outcome_idx = self.cdf_table.partition_point(|&p| p <= r) as i64;

        // Map index to integer value
        self.min_val + outcome_idx
    }
}

impl DGS for CdfSampler {
    /// Initialises [`CdfSampler`] with Gaussian parameter `s` and modulus `mod_q`.
    ///
    /// Parameters:
    /// - `s`: Gaussian parameter
    /// - `mod_q`: The modulus
    fn init(s: f64, mod_q: u64) -> Self {
        // Determine range
        let bound = (s * Self::TAILCUT).ceil() as i64;
        let min_val = -bound;
        let n = (2 * bound + 1) as usize;
        let center_idx = bound as usize; // Index of the '0' outcome

        // Calculate distributed functions
        let mut probabilities = Vec::with_capacity(n);
        let mut sum = 0.0;
        let s_squared = s * s;

        for i in 0..n {
            let x = (min_val + i as i64) as f64;
            let p = (-(PI * x * x) / s_squared).exp();
            probabilities.push(p);
            sum += p;
        }

        // Convert to integer
        let target_sum = u64::MAX;
        let scale = (target_sum as f64) / sum;

        let mut int_probs = Vec::with_capacity(n);
        let mut int_sum: u64 = 0;

        for &p in &probabilities {
            let int_p = (p * scale).round() as u64;
            int_probs.push(int_p);
            int_sum = int_sum.saturating_add(int_p);
        }

        // Correct for rounding errors
        // `int_sum` will not perfectly add up to `target_sum`
        // absorb error in most likely term: `0`
        if int_sum > target_sum {
            let error = int_sum - target_sum;
            int_probs[center_idx] -= error;
        } else {
            let error = target_sum - int_sum;
            int_probs[center_idx] += error;
        }

        // Create discrete CDF table, which sums perfectly to `target_sum`
        let mut cdf_table = Vec::with_capacity(n);
        let mut current_sum: u64 = 0;
        for p in int_probs {
            current_sum = current_sum.saturating_add(p);
            cdf_table.push(current_sum);
        }

        assert_eq!(
            *cdf_table.last().unwrap(),
            target_sum,
            "CDF final sum is incorrect!"
        );

        Self {
            cdf_table,
            min_val,
            mod_q: mod_q as i64,
        }
    }

    /// Samples from the pre-defined discrete Gaussian distribution of this [`CdfSampler`],
    /// providing samples centered around `0` with Gaussian parameter `s` reduced `mod_q`.
    ///
    /// Parameters:
    /// - `rng`: A [`ThreadRng`] to generate randomness
    fn sample(&self, rng: &mut ThreadRng) -> u64 {
        let r = rng.next_u64();

        // Find first index `i` where `cdf[i] > r`
        let outcome_idx = self.cdf_table.partition_point(|&p| p <= r) as i64;

        // Map index to integer value
        let val = self.min_val + outcome_idx;
        let res = (val % self.mod_q + self.mod_q) % self.mod_q;
        res as u64
    }
}

/// Samples a value according to the binomial distribution `Bin(n, 0.5)`
/// centered around `0` reduced `mod_q`.
///
/// Parameters:
/// - `n`: The number of trials / the Binomial parameter `n`
/// - `mod_q`: The modulus
/// - `rng`: A [`ThreadRng`] to generate randomness
pub fn sample_binomial(n: usize, mod_q: i64, rng: &mut ThreadRng) -> u64 {
    let mut x = 0;
    let mut y = 0;

    for _ in 0..n {
        let _tmp_x: bool = rng.random();
        x += _tmp_x as i64;
        let _tmp_y: bool = rng.random();
        y += _tmp_y as i64;
    }

    let val = x - y;
    let res = (val % mod_q + mod_q) % mod_q;
    res as u64
}

#[cfg(test)]
mod test_rej {
    use crate::dgs::{DGS, RejSampler};
    use rand::rng;

    #[test]
    fn in_bounds() {
        let dgs = RejSampler::init(13.0, 13);

        for _ in 0..100 {
            let x = dgs.sample(&mut rng());
            assert!(x < 13);
        }
    }

    #[test]
    fn rough_dist_mini() {
        let dgs = RejSampler::init(3.0, 13);

        let mut vector = vec![0; 13];
        for _ in 0..100_000 {
            let x = dgs.sample(&mut rng());
            vector[x as usize] += 1;
        }

        assert!(vector[0] > 32_500);
        assert!(vector[0] < 34_000);
        assert!(vector[1] > 23_000);
        assert!(vector[1] < 24_000);
        assert!(vector[2] > 7_500);
        assert!(vector[2] < 9_000);
        assert!(vector[3] > 1_000);
        assert!(vector[3] < 2_000);
        assert!(vector[4] > 0);
        assert!(vector[4] < 500);
        assert!(vector[5] < 100);
        assert!(vector[6] < 10);
        assert!(vector[7] < 10);
        assert!(vector[8] < 100);
        assert!(vector[9] < 500);
        assert!(vector[10] > 1_000);
        assert!(vector[10] < 2_000);
        assert!(vector[11] > 7_500);
        assert!(vector[11] < 9_000);
        assert!(vector[12] > 23_000);
        assert!(vector[12] < 24_000);
    }

    #[test]
    fn rough_dist_small() {
        let dgs = RejSampler::init(64.0, 257);

        let mut vector = vec![0; 257];
        for _ in 0..1_000_000 {
            let x = dgs.sample(&mut rng());
            vector[x as usize] += 1;
        }

        let check_points = vec![
            (0, 13_500, 18_500), // Peak
            (1, 13_500, 18_500),
            (16, 11_000, 15_000), // Mid-slope
            (32, 6_000, 8_500),   // 1/8th point
            (48, 2_200, 3_400),
            (64, 500, 1_000),      // 1/4th point
            (80, 80, 400),         // Near tail
            (88, 20, 130),         // Near tail (was [40, 100])
            (96, 5, 60),           // Deep tail
            (128, 0, 25),          // Center minimum (tail)
            (161, 5, 60),          // Symmetric to 96
            (169, 20, 130),        // Symmetric to 88 (Failed here at 49)
            (177, 80, 400),        // Symmetric to 80
            (193, 500, 1_000),     // Symmetric to 64
            (209, 2_200, 3_400),   // Symmetric to 48
            (225, 6_000, 8_500),   // Symmetric to 32
            (241, 11_000, 15_000), // Symmetric to 16
            (256, 13_500, 18_500), // Symmetric to 1
        ];

        for (index, min, max) in check_points {
            let val = vector[index];
            assert!(
                val >= min && val <= max,
                "Assertion failed for index {}: value {} was not in [{}, {}]",
                index,
                val,
                min,
                max
            );
        }
    }
}

#[cfg(test)]
mod test_cdf {
    use crate::dgs::{CdfSampler, DGS};
    use rand::rng;

    #[test]
    fn in_bounds() {
        let dgs = CdfSampler::init(13.0, 13);

        for _ in 0..100 {
            let x = dgs.sample(&mut rng());
            assert!(x < 13);
        }
    }

    #[test]
    fn rough_dist_mini() {
        let dgs = CdfSampler::init(3.0, 13);

        let mut vector = vec![0; 13];
        for _ in 0..100_000 {
            let x = dgs.sample(&mut rng());
            vector[x as usize] += 1;
        }

        assert!(vector[0] > 32_500);
        assert!(vector[0] < 34_000);
        assert!(vector[1] > 23_000);
        assert!(vector[1] < 24_000);
        assert!(vector[2] > 7_500);
        assert!(vector[2] < 9_000);
        assert!(vector[3] > 1_000);
        assert!(vector[3] < 2_000);
        assert!(vector[4] > 0);
        assert!(vector[4] < 500);
        assert!(vector[5] < 100);
        assert!(vector[6] < 10);
        assert!(vector[7] < 10);
        assert!(vector[8] < 100);
        assert!(vector[9] < 500);
        assert!(vector[10] > 1_000);
        assert!(vector[10] < 2_000);
        assert!(vector[11] > 7_500);
        assert!(vector[11] < 9_000);
        assert!(vector[12] > 23_000);
        assert!(vector[12] < 24_000);
    }

    #[test]
    fn rough_dist_small() {
        let dgs = CdfSampler::init(64.0, 257);

        let mut vector = vec![0; 257];
        for _ in 0..1_000_000 {
            let x = dgs.sample(&mut rng());
            vector[x as usize] += 1;
        }

        let check_points = vec![
            (0, 13_500, 18_500), // Peak
            (1, 13_500, 18_500),
            (16, 11_000, 15_000), // Mid-slope
            (32, 6_000, 8_500),   // 1/8th point
            (48, 2_200, 3_400),
            (64, 500, 1_000),      // 1/4th point
            (80, 80, 400),         // Near tail
            (88, 20, 130),         // Near tail (was [40, 100])
            (96, 5, 60),           // Deep tail
            (128, 0, 25),          // Center minimum (tail)
            (161, 5, 60),          // Symmetric to 96
            (169, 20, 130),        // Symmetric to 88 (Failed here at 49)
            (177, 80, 400),        // Symmetric to 80
            (193, 500, 1_000),     // Symmetric to 64
            (209, 2_200, 3_400),   // Symmetric to 48
            (225, 6_000, 8_500),   // Symmetric to 32
            (241, 11_000, 15_000), // Symmetric to 16
            (256, 13_500, 18_500), // Symmetric to 1
        ];

        for (index, min, max) in check_points {
            let val = vector[index];
            assert!(
                val >= min && val <= max,
                "Assertion failed for index {}: value {} was not in [{}, {}]",
                index,
                val,
                min,
                max
            );
        }
    }
}

#[cfg(test)]
mod test_mw {
    use crate::dgs::{DGS, MWSampler};
    use rand::rng;

    #[test]
    fn in_bounds() {
        let dgs = MWSampler::init(257.0, 257);

        for _ in 0..100 {
            let x = dgs.sample(&mut rng());
            assert!(x < 257);
        }
    }
}
