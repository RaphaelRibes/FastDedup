[![Pixi Badge](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/prefix-dev/pixi/main/assets/badge/v0.json)](https://pixi.sh)

# FDedup

FDedup (FastDedup) is a fast and memory-efficient FASTX deduplication tool written in Rust. It utilizes `needletail` for high-performance sequence parsing, `xxhash-rust` for rapid hashing, and `fxhash` for a low-overhead memory cache.

## Features

- **Fast & Memory Efficient**: Uses zero-allocation sequence parsing and a non-cryptographic high-speed hashing cache.
- **Supports Compressed Formats**: Transparently reads and writes both uncompressed and GZIP compressed (`.gz`) FASTQ/FASTA files.
- **Incremental Deduplication**: By default, FDedup can append new sequences to an existing deduplicated output file while pre-loading its existing hashes to prevent any duplicates across runs and helps for crash recovery.
- **Workflow Ready**: Includes a `pixi.toml` file with tasks supporting simulated genome generation, deduplication benching, and FastQC/MultiQC reporting.
- **Profiling Built-in**: Easy memory and execution profiling tasks available via `samply`.

## Requirements

- [Rust](https://rustup.rs/) (>= 1.93)
- [Pixi](https://pixi.sh) (Optional, for running workflows and benchmarks)


## Usage

```bash
fdedup <input_file> [output_file] [--force] [--verbose|-v]
```

- `<input_file>`: Path to the input FASTA/FASTQ/GZ file.
- `[output_file]`: Path to the output file (optional). Defaults to `output.fastq.gz`.
- `--force`: Overwrite the output file if it exists (instead of pre-loading hashes and appending).
- `--verbose` or `-v`: Print processing stats, such as execution time, number of sequences, and duplication rates.
  singularity run fdedup.sif fdedup

### Run it from Cargo
You can run it directly from Cargo:

```bash
cargo run --release -- <input_file> [output_file] [--force] [--verbose|-v]
```

### Run with Pixi

You can also rely on Pixi to run:

```bash
pixi run cargo build --release
pixi run fdedup <input_file> [output_file] [--force] [--verbose|-v]
```

### Run with Singularity / Apptainer
A pre-built Singularity image (`fdedup.sif`) is available for immediate use. You can run the application directly through it:

```bash
singularity run fdedup.sif fdedup <input_file> [output_file] [--force] [--verbose|-v]
```


## To-Do List

- [ ] Support for **Paired-End read deduplication**.
- [ ] Add **Multithreading** to parallelize sequence hashing and processing.
- [ ] Support tracking sequence **abundances** (counts) instead of naive discarding.
- [ ] Add an option for exporting sequences as **FASTA**.
- [ ] Maintain paired qualities correctly for more complex file conversions.
- [ ] Improve error handling.
