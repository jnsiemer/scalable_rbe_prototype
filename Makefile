# Default
.PHONY: all
all: install-os-deps install-rust build test bench
	@echo "====================================================="
	@echo "All done."

# ------------------------------------------------------------------
# ENVIRONMENT SETUP (OS-Dependent)
# ------------------------------------------------------------------

# Install OS-level dependencies (requires Debian and sudo privileges)
.PHONY: install-os-deps
install-os-deps:
	@echo "--- Installing OS Dependencies (only works for Debian-based Linux distributions) ---"
	@echo "Note: You may be prompted for your sudo password."
	sudo apt-get update && sudo apt-get install -y build-essential curl

# Install Rust if not installed
.PHONY: install-rust
install-rust:
	@echo "--- Checking for Rust installation (only works for Unix-like operating systems) ---"
	@if ! command -v cargo >/dev/null 2>&1; then \
		echo "Rust not found. Installing Rust..."; \
		curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y; \
	else \
		echo "Rust is already installed!"; \
	fi


# ------------------------------------------------------------------
# RELEVANT COMMANDS (OS-Independent)
# ------------------------------------------------------------------

# Execute all essential commands
.PHONY: verify
verify: build test bench
	@echo "====================================================="
	@echo "Verification complete! Build, tests, and benchmarks finished."

# Build the project
.PHONY: build
build:
	@echo "--- Building the project ---"
	cargo build --release

# Run the tests
.PHONY: test
test:
	@echo "--- Running unit tests ---"
	cargo test --release -- --test-threads=1

# Run the benchmarks
.PHONY: bench
bench:
	@echo "--- Running benchmarks ---"
	cargo bench RBE
	@echo "--- Benchmarks done ---"
	@echo "Please note comments about execution times in the README."

# Clean the project
.PHONY: clean
clean:
	@echo "--- Cleaning the project ---"
	cargo clean
