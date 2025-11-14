//! Utility functions to (de-)compress and byte (de-)compose matrices.

use std::{fs, path::Path};

/// Ensures that there is no data
pub fn drop_db() {
    let path = "./registry_db";
    if Path::new(path).exists() {
        let _ = fs::remove_dir_all(path);
    }
}

/// Compresses the given `vector` by removing / shifting lower order bits
/// a `prune_by` number of bits to the right for each entry.
///
/// Parameters:
/// - `vector`: The vector with values to compress
/// - `prune_by`: Defines the number of lower order bits to shift / remove from each entry
/// - `mod_q`: The modulus
pub fn compress(vector: &mut [u64], prune_by: usize, mod_q: u64) {
    if prune_by == 0 {
        return;
    }
    for v in vector {
        if *v > mod_q / 2 {
            *v = mod_q - *v;

            *v >>= prune_by - 1;
            if *v & 1 == 1 {
                *v = (*v + 1) % (mod_q >> (prune_by - 1));
            }
            *v >>= 1;
            *v = (mod_q >> prune_by) - *v;
        } else {
            *v >>= prune_by - 1;
            if *v & 1 == 1 {
                *v = (*v + 1) % (mod_q >> (prune_by - 1));
            }
            *v >>= 1;
        }
    }
}

/// Decompresses a compressed `vector` by shifting bits `prune_by` number of bits to the left for each entry.
///
/// Parameters:
/// - `vector`: The vector with values to decompress
/// - `prune_by`: Defines the number of bits to shift left for each entry
/// - `mod_q`: The modulus
pub fn decompress(vector: &mut [u64], prune_by: usize, mod_q: u64) {
    if prune_by == 0 {
        return;
    }
    for v in vector {
        if *v > (mod_q >> prune_by) / 2 {
            *v = (mod_q >> prune_by) - *v;
            *v <<= prune_by;
            *v = mod_q - *v;
        } else {
            *v <<= prune_by;
        }
    }
}

/// Byte-decomposes the [`u64`]-`values`.
///
/// Parameters:
/// - `values`: The values to decompose into bytes
/// - `q`: The largest value that can occur in `values`
pub fn pack(values: &[u64], q: u64) -> Vec<u8> {
    assert!(q > 1);

    let bits = (64 - (q - 1).leading_zeros()) as usize;
    let total_bits = values.len() * bits;
    let total_bytes = total_bits.div_ceil(8);

    let mut out = Vec::with_capacity(total_bytes);
    let mut buffer: u64 = 0;
    let mut buffer_bits = 0;

    for &v in values {
        debug_assert!(v < q);
        buffer |= v << buffer_bits;
        buffer_bits += bits;

        while buffer_bits >= 8 {
            out.push(buffer as u8);
            buffer >>= 8;
            buffer_bits -= 8;
        }
    }

    if buffer_bits > 0 {
        out.push(buffer as u8);
    }

    out
}

/// Composes a vector of [`u64`]-values from bytes.
///
/// Parameters:
/// - `bytes`: The bytes to compose into [`u64`]
/// - `q`: The largest value that occured before decomposing
/// - `len`: Number of entries that was decomposed before
pub fn unpack(bytes: &[u8], q: u64, len: usize) -> Vec<u64> {
    let bits = (64 - (q - 1).leading_zeros()) as usize;
    let mask = (1u64 << bits) - 1;

    let mut values = Vec::with_capacity(len);
    let mut buffer: u64 = 0;
    let mut buffer_bits = 0;

    for &b in bytes {
        buffer |= (b as u64) << buffer_bits;
        buffer_bits += 8;

        while buffer_bits >= bits && values.len() < len {
            values.push(buffer & mask);
            buffer >>= bits;
            buffer_bits -= bits;
        }
    }

    values
}

#[cfg(test)]
mod test_compression {
    use crate::utils::{compress, decompress};
    use rand::{Rng, rng};

    #[test]
    fn round_trip() {
        let q = 257;
        let prune_by = 4;

        let mut vector: Vec<u64> = (0..50).map(|_| rng().random_range(0..q)).collect();
        let cmp = vector.clone();

        compress(&mut vector, prune_by, q);

        decompress(&mut vector, prune_by, q);

        for i in 0..vector.len() {
            assert!(cmp[i] as i64 - vector[i] as i64 >= -2i64.pow(prune_by as u32 - 1));
            assert!(cmp[i] as i64 - vector[i] as i64 <= 2i64.pow(prune_by as u32 - 1));
        }
    }
}
