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
cargo criterion RBE
```
or
```bash
cargo bench RBE
```
To benchmark just a single algorithm, we refer to `./benches/rbe.rs`.

# Parameter Generation

The script for parameter generation can be found in the directory `params`. Please follow the instructions in the [README](params/README.md) to execute the script.

## License

This library is distributed under the [Mozilla Public License Version 2.0](LICENSE).
Permissions of this weak copyleft license are conditioned on making the source code of licensed files and modifications of those files available under the same license (or in certain cases, under one of the GNU licenses). Copyright and license notices must be preserved. Contributors provide an express grant of patent rights. However, a larger work using the licensed work may be distributed under different terms and without source code for files added to the larger work.
