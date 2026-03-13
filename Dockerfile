# Use the official Rust image (meets the >= 1.86.0 requirement)
FROM rust:latest

# Set the working directory inside the container
WORKDIR /usr/src/scalable_rbe_prototype

# Install OS-level dependencies required by the project (like `cc`)
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    && rm -rf /var/lib/apt/lists/*

# Clone the repository
RUN git clone https://github.com/jnsiemer/scalable_rbe_prototype.git .

# Build the library
RUN cargo build --release

# Set the default command when the container launches to an interactive shell,
# allowing you to run benchmarking commands like `cargo bench RBE` manually.
CMD ["/bin/bash"]
