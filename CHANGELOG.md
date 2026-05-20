# Changelog

## v1.3.0 — 2025-05-18

### ✨ New Features

- **Multi-threaded pipeline** (`--threads` / `-j`): Producer/consumer pattern with crossbeam channels for parallel hashing and writing. Auto-detect with `-j 0`.
- **In-process gzip** (flate2+gzippy): Replaced external `pigz` subprocess with in-memory `GzDecoder` for input and `gzippy` multi-threaded compression for output. Up to 3.5× faster at high thread counts vs pigz.
- **Abundance tracking** (`--abundance`): Count occurrences per unique sequence and annotate output headers with `;size=N`.
- **Temp-file acceleration** (`--temp-dir`): Avoid re-decompressing gzip input on the abundance second pass by writing a raw temp file.
- **OOM-safe allocations**: All hash tables use `try_reserve` instead of panicking on out-of-memory conditions.
- **Pre-counting** (`--pre-count`): Scan input first for exact hash table sizing, reducing peak RAM by up to 32%.
- **Paired-end dry-run desync fix**: Corrected desynchronization detection in paired-end `--dry-run` mode.
- **Dry-run without output** (`--dry-run` sans `--output`): `--output` devient optionnel quand `--dry-run` est activé. Pratique pour un quick check avant de lancer le vrai run.

### 🔧 Changes

- `src/writer.rs` extracted from `utils.rs` for cleaner module separation.
- README overhaul with comprehensive CLI docs, abundance/MT/OOM/pre-count flags.
- Added `Alice Baniel` to authors.

### 🐛 Bug Fixes

- Fixed PE dry-run desync detection (task2).
- General error handling improvements.

### ⚙️ Build

- pixi.toml: full build pipeline (conda publish, archive, Apptainer container, Docker).
- pixi-build preview with `RUSTFLAGS="-C link-self-contained=yes"`.