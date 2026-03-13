// Copyright 2026 Jan Niklas Siemer
//
// This file is part of scalable_rbe.
//
// scalable_rbe is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! Contains benchmarks for the Registration-Based Encryption scheme.

use criterion::{Criterion, criterion_group};
use rand::{RngExt as _, distr::Uniform, rng};
use scalable_rbe::{
    mat::Mat,
    ntt::MatNTT,
    rbe::{B, D, K, M_B, N, Q, RBE},
    utils::drop_db,
};

/// Benchmark [RBE::setup].
///
/// This benchmark can be run with for example:
/// - `cargo criterion Single\ Setup`
/// - `cargo bench --bench benchmarks Single\ Setup`
/// - `cargo flamegraph --bench benchmarks -- --bench Single\ Setup`
fn bench_setup(c: &mut Criterion) {
    let rbe = RBE::init(false);

    c.bench_function("Single Setup", |b| b.iter(|| rbe.setup()));

    drop_db();
}

/// Benchmark [RBE::keygen].
///
/// This benchmark can be run with for example:
/// - `cargo criterion Single\ KeyGen`
/// - `cargo bench --bench benchmarks Single\ KeyGen`
/// - `cargo flamegraph --bench benchmarks -- --bench Single\ KeyGen`
fn bench_keygen(c: &mut Criterion) {
    let rbe = RBE::init(false);
    let (pp, _st, mut _reg) = rbe.setup();

    c.bench_function("Single KeyGen", |b| b.iter(|| rbe.keygen(&pp)));

    drop_db();
}

/// Benchmark [RBE::update].
///
/// This benchmark can be run with for example:
/// - `cargo criterion Single\ Update`
/// - `cargo bench --bench benchmarks Single\ Update`
/// - `cargo flamegraph --bench benchmarks -- --bench Single\ Update`
fn bench_update(c: &mut Criterion) {
    let rbe = RBE::init(false);
    let (pp, _st, mut reg) = rbe.setup();
    let (pk, _sk) = rbe.keygen(&pp);

    c.bench_function("Single Update", |b| {
        b.iter(|| rbe.update(&mut reg, &pp, 0, &"0".repeat(162), &pk))
    });

    drop_db();
}

/// Benchmark [RBE::enc] with a discrete Gaussian distribution as error distribution.
///
/// This benchmark can be run with for example:
/// - `cargo criterion Single\ Enc`
/// - `cargo bench --bench benchmarks Single\ Enc`
/// - `cargo flamegraph --bench benchmarks -- --bench Single\ Enc`
fn bench_enc(c: &mut Criterion) {
    let id = "0".repeat(162);

    let rbe = RBE::init(false);
    let (pp, _st, mut reg) = rbe.setup();
    let (pk, _sk) = rbe.keygen(&pp);
    let st = rbe.update(&mut reg, &pp, 0, &id, &pk);
    let msg_dist = Uniform::new(0, 2).unwrap();
    let vec_m_t = Mat::sample_small(D, 1, B, &msg_dist, 0, Q, &mut rng());

    c.bench_function("Single Enc", |b| {
        b.iter(|| rbe.enc(&pp, &st, &id, vec_m_t.clone()))
    });

    drop_db();
}

/// Benchmark [RBE::enc_wo_randomness].
///
/// This benchmark can be run with for example:
/// - `cargo criterion Single\ Enc\ without\ randomness`
/// - `cargo bench --bench benchmarks Single\ Enc\ without\ randomness`
/// - `cargo flamegraph --bench benchmarks -- --bench Single\ Enc\ without\ randomness`
fn bench_enc_wo_randomness(c: &mut Criterion) {
    let id = "0".repeat(162);
    let l = id.len();

    let rbe = RBE::init(false);
    let (pp, _st, mut reg) = rbe.setup();
    let (pk, _sk) = rbe.keygen(&pp);
    let st = rbe.update(&mut reg, &pp, 0, &id, &pk);
    let msg_dist = Uniform::new(0, 2).unwrap();
    let vec_m_t = Mat::sample_small(D, 1, B, &msg_dist, 0, Q, &mut rng());

    // r_j ← χ_s^n ∀j ∈ {0,1,...,ℓ}
    let vecs_r_j_t: Vec<MatNTT> = (0..=l)
        .into_iter()
        .map(|_| {
            let vec_r_j_t = Mat::sample_disc_gauss(D, 1, N, &rbe.dgs_sigma, &mut rng());
            MatNTT::ntt(vec_r_j_t, &rbe.ntt_plan)
        })
        .collect();

    // e_j ← χ_e^{k * m_A} ∀j ∈ [ℓ]
    let vecs_e_j_t: Vec<Mat> = (0..l)
        .into_iter()
        .map(|_| Mat::sample_disc_gauss(D, 1, K * rbe.m_a, &rbe.dgs_sigma, &mut rng()))
        .collect();

    // e_ℓ ← χ_e^{m_B}
    let vec_e_l_t = Mat::sample_disc_gauss(D, 1, M_B, &rbe.dgs_sigma, &mut rng());

    // e¯ ←$ χ¯^T
    // error term s_e is a guess here, not determined according to the document
    let vec_e_tilde_t = Mat::sample_disc_gauss(D, 1, st.nr_columns, &rbe.dgs_e_tilde, &mut rng());

    // let (vecs_r_j_t, vecs_e_j_t, vec_e_l_t, vec_e_tilde_t) = rbe.gen_enc_randomness(162, 1);

    c.bench_function("Single Enc without randomness", |b| {
        b.iter(|| {
            rbe.enc_wo_randomness(
                &pp,
                &st,
                &id,
                vec_m_t.clone(),
                &vecs_r_j_t,
                &vecs_e_j_t,
                &vec_e_l_t,
                &vec_e_tilde_t,
            )
        })
    });

    drop_db();
}

/// Benchmark [RBE::witness_gen].
///
/// This benchmark can be run with for example:
/// - `cargo criterion Single\ Witness\ Gen`
/// - `cargo bench --bench benchmarks Single\ Witness\ Gen`
/// - `cargo flamegraph --bench benchmarks -- --bench Single\ Witness\ Gen`
fn bench_witness_gen(c: &mut Criterion) {
    let id = "0".repeat(162);

    let rbe = RBE::init(false);
    let (pp, _st, mut reg) = rbe.setup();
    let (pk, _sk) = rbe.keygen(&pp);
    let _st = rbe.update(&mut reg, &pp, 0, &id, &pk);

    c.bench_function("Single Witness Gen", |b| {
        b.iter(|| rbe.witness_gen(&reg, 0, &id))
    });

    drop_db();
}

/// Benchmark [RBE::dec].
///
/// This benchmark can be run with for example:
/// - `cargo criterion Single\ Dec`
/// - `cargo bench --bench benchmarks Single\ Dec`
/// - `cargo flamegraph --bench benchmarks -- --bench Single\ Dec`
fn bench_dec(c: &mut Criterion) {
    let id = "0".repeat(162);

    let rbe = RBE::init(false);
    let (pp, _st, mut reg) = rbe.setup();
    let (pk, sk) = rbe.keygen(&pp);
    let st = rbe.update(&mut reg, &pp, 0, &id, &pk);
    let msg_dist = Uniform::new(0, 2).unwrap();
    let vec_m_t = Mat::sample_small(D, 1, B, &msg_dist, 0, Q, &mut rng());

    let ct = rbe.enc(&pp, &st, &id, vec_m_t);
    let wit = rbe.witness_gen(&reg, 0, &id).unwrap();

    c.bench_function("Single Dec", |b| b.iter(|| rbe.dec(&sk, 0, &wit, &ct)));

    drop_db();
}

/// Benchmark [RBE::setup], [RBE::keygen], [RBE::update], [RBE::enc], [RBE::witness_gen], and [RBE::dec].
///
/// This benchmark can be run with for example:
/// - `cargo criterion RBE`
/// - `cargo bench --bench benchmarks RBE`
fn bench_all(c: &mut Criterion) {
    let rbe = RBE::init(false);

    // Benchmark setup
    c.bench_function("RBE Setup", |b| b.iter(|| rbe.setup()));

    // Setup to continue with
    let (pp, _st, mut reg) = rbe.setup();

    // Benchmark keygen
    c.bench_function("RBE KeyGen", |b| b.iter(|| rbe.keygen(&pp)));

    // Generate (pk, sk)-pair that we reuse from now on
    let (pk, sk) = rbe.keygen(&pp);
    let mut ids = vec![];
    let mut rng = rng();

    // Benchmark update
    c.bench_function("RBE Update", |b| {
        b.iter_batched(
            || {
                // This setup phase is not timed in the benchmark.
                let id: String = (0..162).map(|_| rng.random_range('0'..='2')).collect();
                ids.push(id.clone());
                let t = (ids.len() as f64).log2().floor() as usize;
                (id, t)
            },
            |(id, t)| rbe.update(&mut reg, &pp, t, &id, &pk),
            criterion::BatchSize::SmallInput,
        )
    });

    // Get current state
    let st = rbe.update(&mut reg, &pp, 0, &ids[0], &pk);

    // Generate small random message
    let msg_dist = Uniform::new(0, 2).unwrap();
    let vec_m_t = Mat::sample_small(D, 1, B, &msg_dist, 0, Q, &mut rng);

    // Benchmark enc
    c.bench_function("RBE Enc", |b| {
        b.iter_batched(
            || vec_m_t.clone(),
            |vec_m_t| rbe.enc(&pp, &st, &ids[0], vec_m_t),
            criterion::BatchSize::SmallInput,
        )
    });

    // Get ciphertext to decrypt later
    let ct = rbe.enc(&pp, &st, &ids[0], vec_m_t);

    // Benchmark witness_gen
    c.bench_function("RBE Witness Gen", |b| {
        b.iter(|| rbe.witness_gen(&reg, 0, &ids[0]))
    });

    // Get witness
    let wit = rbe.witness_gen(&reg, 0, &ids[0]).unwrap();

    // Benchmark dec
    c.bench_function("RBE Dec", |b| b.iter(|| rbe.dec(&sk, 0, &wit, &ct)));

    // Clear registry database
    drop_db();
}

criterion_group! {
    name = benches;
    config = Criterion::default();
    targets = bench_setup,
        bench_keygen,
        bench_update,
        bench_enc,
        bench_enc_wo_randomness,
        bench_witness_gen,
        bench_dec,
        bench_all,
}
