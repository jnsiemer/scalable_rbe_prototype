// Copyright 2026 Jan Niklas Siemer
//
// This file is part of scalable_rbe.
//
// scalable_rbe is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! # Prototype of 'Scalable Registration-Based Encryption from Lattices'
//! This library contains a prototype of a Registration-Based Encryption (RBE) scheme.
//! The post-quantum security of the implemented construction is based
//! on the standard Learning with Errors assumption.
//! The analysis enables several tweaks to significantly reduce ciphertext
//! sizes in practical deployments. These can be achieved via small extensions of this prototype.
//!
//! The only operation that is not implemented by this prototype is the merge of two
//! SIS-trees. Instead of adding public keys of identities into trees blindly, we
//! insert them with the expected layout in our benchmarks.
//!
//! The layout of the library should provide an API that allows for quick and easy experiments.
//! Thus, several fields and algorithms are `public`, which should not be public in a real application.
//!
//! **WARNING:** This is a prototype. It should not be deployed in real applications.
//! More specifically, it does not provide constant-time security.
//!
//! ## External Libraries / Crates
//! The implementation of this library makes use of Zama's [`concrete_ntt`] crate,
//! which supports hardware accelerated NTT transforms and multiplication if `AVX2` or `AVX512`
//! is supported (using the nightly toolchain for the latter).
//! Furthermore, it uses the crates [`rand`] and [`rand_distr`] to generate 'cryptographically
//! secure' randomness.
//! We use the key-value database [`sled`] to store vectors in the registry and [`rayon`] to
//! enable multi-threading in [`rbe::RBE::enc`].
//! All benchmarks can be found in the directory `benches` and we rely on [`criterion`] to
//! benchmark the implementation.

pub mod dgs;
pub mod mat;
pub mod ntt;
pub mod rbe;
pub mod utils;
