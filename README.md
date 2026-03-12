# Prototype of 'Scalable Registration-Based Encryption from Lattices'
[<img alt="github" src="https://img.shields.io/badge/github-github?style=for-the-badge&logo=github&color=8da0cb" height="20">](https://github.com/jnsiemer/scalable_rbe_prototype)
[<img alt="crates.io" src="https://img.shields.io/badge/crates-cratesio?style=for-the-badge&logo=rust&color=fc8d62" height="20">](https://crates.io/crates/scalable-rbe)
[<img alt="crates.io" src="https://img.shields.io/badge/zenodo-zenodo?style=for-the-badge&logo=zenodo&label=DOI&logoColor=000&color=ffd92f" height="20">](https://doi.org/10.5281/zenodo.18989489)
[<img alt="docs.rs" src="https://img.shields.io/badge/docs-docs?style=for-the-badge&logo=docs.rs&color=66c2a5" height="20">](https://docs.rs/scalable-rbe)
[<img alt="build" src="https://img.shields.io/github/actions/workflow/status/jnsiemer/scalable_rbe_prototype/pipeline.yml?style=for-the-badge" height="20">](https://github.com/jnsiemer/scalable_rbe_prototype/actions/workflows/pipeline.yml)
[<img alt="license" src="https://img.shields.io/badge/License-MPL_2.0-blue.svg?style=for-the-badge" height="20">](LICENSE)

## Dependencies
- [`Rust`](https://rust-lang.org/) version >= 1.86.0 (optional: `nightly` toolchain to enable AVX-512 hardware acceleration if available by passing `--features=nightly` to `cargo bench`)
- Installation of several standard libraries such as `cc`, which come out-of-the-box with several Linux distributions but can be installed using `sudo apt install build-essential` 

To check the version and toolchain, use `rustc --version`. To update `rustup update`. To install the `nightly` toolchain, `rustup toolchain install nightly`.

# Relevant Commands
To build the library, run
```bash
cargo build --release
```

The documentation is available at [docs.rs](https://docs.rs/scalable-rbe) or you can generate it yourself by executing
```bash
cargo doc --open
```
If you're using WSL and your browser does not open, retry after installing `wslu`.

To run all unit-tests, which also ensure correctness of the parameters and construction, execute the following line.
```bash
cargo test --release -- --test-threads=1 # single threaded to prevent simoultaneous writes to the database from different tests
```

To benchmark all RBE algorithms, you can use
```bash
cargo bench RBE
```
If you have [`criterion`](https://crates.io/crates/criterion) installed, you can also execute
```bash
cargo criterion RBE
```
To benchmark just a single algorithm, we refer to `./benches/rbe.rs`.

# Parameter Generation

The script for parameter generation can be found in the directory `params`. Please follow the instructions in the [README](params/README.md) to execute the script.

## License

This library is distributed under the [Mozilla Public License Version 2.0](LICENSE).
Permissions of this weak copyleft license are conditioned on making the source code of licensed files and modifications of those files available under the same license (or in certain cases, under one of the GNU licenses). Copyright and license notices must be preserved. Contributors provide an express grant of patent rights. However, a larger work using the licensed work may be distributed under different terms and without source code for files added to the larger work.
