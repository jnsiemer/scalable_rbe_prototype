# Prototype of 'Scalable Registration-Based Encryption from Lattices'

## Dependencies
- [`Rust`](https://rust-lang.org/) version >= 1.67.0 with the `nightly` toolchain (for AV512 hardware acceleration)

To check the version and toolchain, use `rustc --version`. To update `rustup update`. To install the `nightly` toolchain, `rustup toolchain install nightly`.

# Relevant Commands
To build the library, run
```bash
cargo build --release
```

To generate the documentation, execute
```bash
cargo doc --open
```

To run all unit-tests, which also ensure correctness of the parameters and construction, execute the following line.
```bash
cargo test --release -- --test-threads=1 # single threaded to prevent simoultaneous writes to the database from different tests
```

To benchmark all RBE algorithms, use either
```bash
cargo criterion RBE:
```
or
```bash
cargo bench RBE:
```
To benchmark just a single algorithm, we refer to `./benches/rbe.rs`.

# Parameter Generation
The script for parameter generation can be found in the directory `params`.
