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

## Building on this project
To integrate this library in your project, you can create a dependency via [crates.io](https://crates.io/crates/scalable-rbe) the command
```bash
cargo add scalable-rbe
```
Then, you can use the exposed API `RBE::init(false)`, `rbe.setup()`, `rbe.keygen(&pp)`, etc. to utilise Registration-Based Encryption in your protocol.
> **WARNING:** This is a prototype and should not be used in a production envirnoment. This implementation does _not_ provide any side-channel security.

If you need to use different parameters, please clone/fork this project and create a [dependency](https://doc.rust-lang.org/cargo/reference/specifying-dependencies.html) based on your GitHub repository or local directory. Please ensure that any diverging set of parameters meets the correctness and security requirements specified in our paper. You can use our [parameter selection script](params/param_selection.ipynb) to find valid sets of parameters.

## Executing Tests and Benchmarks

If you want to verify the results of this project, i.e. build the project, run its unit tests, build its documentation, or run the benchmarks, there are three options.
1. Local installation of [dependencies](#dependencies) to enable the local execution of all [relevant commands](#relevant-commands).
2. If you are using a Debian-based operating system, you may utilise [`make`](#makefile) to install all dependencies and run all commands automatically.
3. Follow the [instructions](#docker) to setup a docker container and execute all commands in the container.

### Makefile
If you are running a Debian-based operating system, you may use `make all` to ensure all [dependencies](#dependencies) are properly installed.
Otherwise, you can still use `make verify` to build, test, and run benchmarks for this project and we refer to the [Docker image](#docker) as an alternative to create a container with preinstalled dependencies.

### Docker
We assume that [Docker Desktop](https://www.docker.com/products/docker-desktop/) is setup on your machine.
Then, you can build the appropriate docker image by executing
```bash
docker build -t rbe-env .
```
Once the build is finished, your docker image is saved as `rbe-env` and we can enter the image by running the following command.
```bash
docker run -it rbe-env
```
Now, you can run `make verify` to build, test, and run benchmarks for this project.

Alternatively, you can interactively run any `cargo` commands inside your Docker container. Use `exit` to leave the container.

### Relevant Commands
Ensure that all [dependencies](#dependencies) are setup as required. Then, the following yields a list of valid commands to build, test, and benchmark this project.

To build the library, run
```bash
cargo build --release
```

The documentation is available at [docs.rs](https://docs.rs/scalable-rbe) or you can generate it yourself by executing
```bash
cargo doc --open
```
If you're using [WSL](https://github.com/microsoft/WSL) and your browser does not open, retry after installing `wslu`.

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

#### Note on Benchmarks
The benchmarks in the paper were conducted on an Intel Core Ultra 7 165U 4.9 GHz CPU with 32 GB of RAM using AVX2 / SIMD acceleration and 14 threads.
The number of threads impacts the execution time of encryption `Enc` significantly (as it is the only multi-threaded algorithm).
Our experimental results show that disabling AVX2 instructions may slow down the execution time of some algorithms up to a factor of 2.
We did not conduct any experiments with AVX-512 acceleration.

---

# Parameter Generation

The script for parameter generation can be found in the directory `params`. Please follow the instructions in the [README](params/README.md) to execute the script.

## License

This library is distributed under the [Mozilla Public License Version 2.0](LICENSE).
Permissions of this weak copyleft license are conditioned on making the source code of licensed files and modifications of those files available under the same license (or in certain cases, under one of the GNU licenses). Copyright and license notices must be preserved. Contributors provide an express grant of patent rights. However, a larger work using the licensed work may be distributed under different terms and without source code for files added to the larger work.
