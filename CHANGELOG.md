# Changelog

All notable changes to this project will be documented in this file.

## [1.2.0] - Unreleased

### Added
- **Library API (`lib.rs`)**: FDedup is now completely modularized and can be used as a standalone Rust library!
  - Introduced `UniqueFastxStream` for high-performance, deterministic streaming deduplication using zero-copy callbacks (`for_each_unique`).
  - Added programmatic capacity hints (`with_capacity`) to pre-allocate memory and eliminate rehashing.
  - Exposed 64-bit (`u64`) and 128-bit (`u128`) explicit hashing traits directly to library consumers avoiding CLI overhead.
  - Generates comprehensive Library Usage documentation within `README.md`.
- **Cargo Modularization**: Made the `clap` CLI parser an optional dependency (`cli` feature). This dramatically reduces the dependency tree footprint when compiling FastDedup strictly as a library crate.
- Added the `--read-length` (`-P`) argument to CLI to enable tuning of internal I/O write buffer sizes and sequence capacity estimates based on the expected read length in base pairs (defaults to 150bp).
- Implemented `compute_write_buffer_size` and `compute_bytes_per_read` methods to optimally constrain memory buffer operations contextually based on read formats.
- Added exhaustive testing coverage for buffer estimations, CLI options scaling, and trait-based hash determinism (`test_hash_sequence_deterministic`, `test_hash_pair_order_dependent`, etc).

### Changed
- Refactored `HashMode` matching in CLI into standard `.from()` Traits mapping directly to the now-public `HashType` enum structure.
- Translated all test functions, assertions, and inline documentation in `hasher.rs` and `utils.rs` from French to English for broader developer accessibility (e.g., `test_birthday_64bit_a_2_puissance_32` -> `test_birthday_64bit_at_2_pow_32`).
- Streamlined FastA serializing logic via standardizing around a new reusable `write_fasta_record` utility method.

### Fixed
- **Buffer Safety**: Enforced output buffer flushing (`buf_writer.flush()`) explicitly at the conclusion of both single-end and paired-end write tasks to completely flush I/O streams safely before exit.
- Enforced stricter CLI argument validations using Clap: `--output-r2` (`-p`) is explicitly required by `--input-r2` (`-2`) at the parser level, preventing edge-case crashes.
