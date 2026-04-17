use anyhow::{bail, Context, Result};
use needletail::parse_fastx_file;
use std::io::Write;
use std::path::Path;
use std::{fs, io};
use std::fs::{File, OpenOptions};
use flate2::write::GzEncoder;
use flate2::Compression;
use crate::hasher::{SequenceHasher, HashType, HashVerifier};

/// Maximum number of entries to pre-allocate in the hash set.
/// Prevents multi-GB allocations from inaccurate file-size estimates.
const MAX_PREALLOC_ENTRIES: usize = 500_000_000; // 500M

/// Supported output formats for sequences.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum OutputFormat {
    Fasta,
    Fastq,
}

impl OutputFormat {
    /// Determines the output format based on the file extension.
    pub fn from_extension(path: &Path) -> Self {
        let path_str = path.to_string_lossy().to_lowercase();

        let without_gz = path_str.strip_suffix(".gz").unwrap_or(&path_str);

        if without_gz.ends_with(".fasta") || without_gz.ends_with(".fa") || without_gz.ends_with(".fna") {
            OutputFormat::Fasta
        } else {
            OutputFormat::Fastq
        }
    }

    /// Checks if the file path indicates a gzipped file.
    pub fn is_gz(path: &Path) -> bool {
        path.to_string_lossy().to_lowercase().ends_with(".gz")
    }
}

pub fn birthday_problem_square_approximation(x: usize, n: &HashType) -> f64 {
    (x as f64).powi(2) / 2.0_f64.powi(n.to_num() as i32 + 1)
}

/// Returns the appropriate `HashType` given the dataset size and collision probability threshold.
pub fn get_hash_method(size: usize, threshold: f64) -> HashType {
    if (2f64 * 2.0f64.powi(64) * threshold).sqrt() < size as f64 {
        HashType::XXH3_128
    } else {
        HashType::XXH3_64
    }
}

/// Computes the write buffer size from the expected read length.
///
/// A FASTQ record is roughly:
///   - 1 ID line:      ~40 bytes (conservative estimate)
///   - 1 sequence:     read_length bytes
///   - 1 '+' line:     2 bytes
///   - 1 quality line: read_length bytes
///
/// Total ≈ `2 * read_length + 42` bytes.
///
/// We target holding ~512 records per flush as a heuristic — large enough
/// to amortize syscall and gzip overhead, small enough to avoid memory pressure
/// when many instances run in parallel on an HPC node.
/// The result is clamped between 64 KB (floor, never regress below OS page granularity)
/// and 64 MB (ceiling, avoids OOM on ultra-long reads at high parallelism).
pub fn compute_write_buffer_size(read_length: usize) -> usize {
    let bytes_per_record = 2 * read_length + 42;
    let target = bytes_per_record * 512;
    target.clamp(64 * 1024, 64 * 1024 * 1024)
}

/// Computes the bytes-per-read divisor for capacity pre-allocation.
///
/// Same record size formula as above, with a compression factor of ~4x for gzip.
/// Clamped to a minimum of 1 to prevent division by zero.
pub fn compute_bytes_per_read(read_length: usize, is_gz: bool) -> usize {
    let plain = (2 * read_length + 42).max(1);
    if is_gz {
        (plain / 4).max(1) // ~3.5–4x compression on DNA is typical
    } else {
        plain
    }
}

/// Prepares a writer for the output file based on the format and whether it should
/// be forced (overwritten) or appended.
///
/// # Errors
///
/// Returns an error if:
/// - Append mode is requested on a `.gz` file (multi-member gzip is unreliable
///   for downstream bioinformatics tools).
/// - The file cannot be opened or created.
pub fn prepare_writer(path: &Path, force: bool, compression_level: u32) -> Result<(Box<dyn Write>, OutputFormat)> {
    let format = OutputFormat::from_extension(path);
    let is_gz = OutputFormat::is_gz(path);

    let file = if path.exists() && !force {
        if is_gz {
            bail!(
                "Cannot append to gzipped output {:?}. Multi-member gzip is unreliable for \
                 downstream tools. Use --force to overwrite, or use an uncompressed output \
                 format for incremental/resumable runs.",
                path
            );
        }
        OpenOptions::new().append(true).open(path)
            .with_context(|| format!("Failed to open file in append mode: {:?}", path))?
    } else {
        File::create(path)
            .with_context(|| format!("Failed to create file: {:?}", path))?
    };

    let writer: Box<dyn Write> = if is_gz {
        Box::new(GzEncoder::new(file, Compression::new(compression_level)))
    } else {
        Box::new(file)
    };

    Ok((writer, format))
}

/// Estimates the number of sequences in a file to pre-allocate memory.
/// The result is capped at `MAX_PREALLOC_ENTRIES` to avoid over-allocation.
pub fn estimate_sequence_capacity<P: AsRef<Path>>(
    path: P,
    read_length: usize,
) -> Result<usize> {
    let path_ref = path.as_ref();
    if !path_ref.exists() {
        return Ok(0);
    }

    let file_size_bytes = fs::metadata(path_ref)?.len() as usize;
    let is_gz = OutputFormat::is_gz(path_ref);
    let divisor = compute_bytes_per_read(read_length, is_gz);

    Ok((file_size_bytes / divisor).clamp(1, MAX_PREALLOC_ENTRIES))
}

/// Preloads existing hashes into memory from an existing output file (Single-End).
pub fn preload_existing_hashes<T: SequenceHasher>(
    path: &str,
    verifier: &mut HashVerifier<T>,
    verbose: bool,
) -> Result<(usize, u64)> {
    if !Path::new(path).exists() {
        return Ok((0, 0));
    }

    if verbose {
        println!("Preloading sequences from existing output...");
    }

    let mut reader = parse_fastx_file(path)
        .context("Error opening preload file")?;
    let mut count = 0;
    let mut valid_bytes = ByteCounter(0);

    while let Some(record_result) = reader.next() {
        match record_result {
            Ok(record) => {
                let hash = T::hash_sequence(&record.seq());
                verifier.verify(hash);

                count += 1;
                let _ = record.write(&mut valid_bytes, None);
            }
            Err(e) => {
                if verbose {
                    eprintln!("Incomplete sequence detected at end of file ({}).", e);
                    eprintln!("Calculating fail-safe truncation point...");
                }
                break;
            }
        }
    }

    Ok((count, valid_bytes.0))
}

/// Preloads existing hashes into memory from existing output files (Paired-End)
/// and computes valid byte sizes to synchronize truncation.
pub fn preload_existing_paired_hashes<T: SequenceHasher>(
    path_r1: &str,
    path_r2: &str,
    verifier: &mut HashVerifier<T>,
    verbose: bool,
) -> Result<(usize, u64, u64)> {
    let path_r1_ref = Path::new(path_r1);
    let path_r2_ref = Path::new(path_r2);

    if !path_r1_ref.exists() || !path_r2_ref.exists() {
        return Ok((0, 0, 0));
    }

    if verbose {
        println!("Preloading and synchronizing pairs from existing outputs...");
    }

    let mut reader_r1 = parse_fastx_file(path_r1).context("Error opening R1 preload file")?;
    let mut reader_r2 = parse_fastx_file(path_r2).context("Error opening R2 preload file")?;

    let mut count = 0;
    let mut valid_bytes_r1 = ByteCounter(0);
    let mut valid_bytes_r2 = ByteCounter(0);

    while let (Some(record_r1_res), Some(record_r2_res)) = (reader_r1.next(), reader_r2.next()) {
        match (record_r1_res, record_r2_res) {
            (Ok(record_r1), Ok(record_r2)) => {
                let combined_hash = T::hash_pair(&record_r1.seq(), &record_r2.seq());
                verifier.verify(combined_hash);

                count += 1;
                let _ = record_r1.write(&mut valid_bytes_r1, None);
                let _ = record_r2.write(&mut valid_bytes_r2, None);
            }
            _ => {
                if verbose {
                    eprintln!("Incomplete sequence or desynchronization detected at end of file.");
                    eprintln!("Calculating fail-safe truncation points for R1 and R2...");
                }
                break;
            }
        }
    }

    Ok((count, valid_bytes_r1.0, valid_bytes_r2.0))
}

/// A no-op writer that counts bytes written, used to compute valid byte offsets
/// without actually performing I/O.
pub(crate) struct ByteCounter(pub u64);

impl Write for ByteCounter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.0 += buf.len() as u64;
        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

/// Writes a single FASTA record using raw byte.
#[inline]
pub fn write_fasta_record(writer: &mut impl Write, id: &[u8], seq: &[u8]) -> Result<()> {
    writer.write_all(b">").context("Error writing FASTA '>' prefix")?;
    writer.write_all(id).context("Error writing FASTA ID")?;
    writer.write_all(b"\n").context("Error writing FASTA ID newline")?;
    writer.write_all(seq).context("Error writing FASTA sequence")?;
    writer.write_all(b"\n").context("Error writing FASTA sequence newline")?;
    Ok(())
}

/// Validates that the output format is compatible with the input data.
/// FASTA input cannot be written as FASTQ (no quality scores).
pub fn validate_format_compatibility(has_quality: bool, output_format: OutputFormat) -> Result<()> {
    if !has_quality && output_format == OutputFormat::Fastq {
        bail!(
            "FASTA → FASTQ conversion not supported: the input file does not contain \
             quality scores. Use a .fasta/.fa/.fna extension for output."
        );
    }
    Ok(())
}

/// Extracts the base read ID, stripping both the description (anything after
/// the first whitespace) and the `/1` or `/2` mate suffix used by pre-CASAVA-1.8
/// Illumina headers and some SRA `fastq-dump --split-files` outputs.
#[inline]
pub fn base_read_id(id: &[u8]) -> &[u8] {
    // Strip description: "ID desc" → "ID"
    let head = id.split(|&b| b == b' ' || b == b'\t').next().unwrap_or(id);
    // Strip /1 or /2 mate suffix if present
    if head.len() >= 2 {
        let tail = &head[head.len() - 2..];
        if tail == b"/1" || tail == b"/2" {
            return &head[..head.len() - 2];
        }
    }
    head
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hasher::HashType;

    // --- birthday_problem_square_approximation ---

    #[test]
    fn test_birthday_64bit_at_2_pow_32() {
        let n = (u32::MAX as usize) + 1;
        let prob = birthday_problem_square_approximation(n, &HashType::XXH3_64);
        assert!(
            (prob - 0.5).abs() < 1e-9,
            "Expected probability 0.5, got {:.12}",
            prob
        );
    }

    #[test]
    fn test_birthday_128bit_at_2_pow_32() {
        let n = (u32::MAX as usize) + 1;
        let prob = birthday_problem_square_approximation(n, &HashType::XXH3_128);
        let expected = (n as f64).powi(2) / 2.0_f64.powi(129);
        assert!(
            (prob - expected).abs() < 1e-30,
            "Expected probability {:.3e}, got {:.3e}",
            expected,
            prob
        );
    }

    #[test]
    fn test_birthday_zero_sequences() {
        let prob = birthday_problem_square_approximation(0, &HashType::XXH3_64);
        assert_eq!(prob, 0.0);
    }

    #[test]
    fn test_birthday_one_sequence() {
        let prob = birthday_problem_square_approximation(1, &HashType::XXH3_64);
        assert!(prob > 0.0 && prob < 1e-15);
    }

    // --- get_hash_method ---

    #[test]
    fn test_get_hash_method_large_dataset_selects_128() {
        let large_size = 200_000_000_usize;
        let method = get_hash_method(large_size, 0.001);
        assert!(
            matches!(method, HashType::XXH3_128),
            "A dataset of {} sequences should use XXH3_128", large_size
        );
    }

    #[test]
    fn test_get_hash_method_small_dataset_selects_64() {
        let small_size = 1_000_usize;
        let method = get_hash_method(small_size, 0.001);
        assert!(
            matches!(method, HashType::XXH3_64),
            "A dataset of {} sequences should use XXH3_64", small_size
        );
    }

    #[test]
    fn test_get_hash_method_threshold_zero_forces_128() {
        let method = get_hash_method(1, 0.0);
        assert!(
            matches!(method, HashType::XXH3_128),
            "A threshold of 0 should always select XXH3_128"
        );
    }

    #[test]
    fn test_get_hash_method_threshold_one_forces_64() {
        let method = get_hash_method(1_000_000, 1.0);
        assert!(
            matches!(method, HashType::XXH3_64),
            "A threshold of 1.0 should select XXH3_64 for a normal-sized dataset"
        );
    }

    #[test]
    fn test_get_hash_method_tipping_point() {
        let boundary = 4_294_967_296_usize;
        let just_below = get_hash_method(boundary - 1, 0.5);
        let just_above = get_hash_method(boundary + 1, 0.5);
        assert!(
            matches!(just_below, HashType::XXH3_64),
            "Just below the boundary should be XXH3_64"
        );
        assert!(
            matches!(just_above, HashType::XXH3_128),
            "Just above the boundary should be XXH3_128"
        );
    }

    // --- compute_write_buffer_size ---

    #[test]
    fn test_buffer_size_illumina_150bp() {
        let buf = compute_write_buffer_size(150);
        assert_eq!(buf, 175_104);
        assert!(buf >= 64 * 1024);
        assert!(buf <= 64 * 1024 * 1024);
    }

    #[test]
    fn test_buffer_size_floor_clamping() {
        let buf = compute_write_buffer_size(1);
        assert_eq!(buf, 64 * 1024);
    }

    #[test]
    fn test_buffer_size_ceiling_clamping() {
        let buf = compute_write_buffer_size(100_000);
        assert_eq!(buf, 64 * 1024 * 1024);
    }

    // --- compute_bytes_per_read ---

    #[test]
    fn test_bytes_per_read_plain_150bp() {
        assert_eq!(compute_bytes_per_read(150, false), 342);
    }

    #[test]
    fn test_bytes_per_read_gz_150bp() {
        assert_eq!(compute_bytes_per_read(150, true), 85);
    }

    #[test]
    fn test_bytes_per_read_zero_length_no_panic() {
        assert!(compute_bytes_per_read(0, false) >= 1);
        assert!(compute_bytes_per_read(0, true) >= 1);
    }

    // --- OutputFormat ---

    #[test]
    fn test_output_format_fastq_default() {
        assert_eq!(OutputFormat::from_extension(Path::new("out.fastq")), OutputFormat::Fastq);
        assert_eq!(OutputFormat::from_extension(Path::new("out.fq.gz")), OutputFormat::Fastq);
        assert_eq!(OutputFormat::from_extension(Path::new("out.unknown")), OutputFormat::Fastq);
    }

    #[test]
    fn test_output_format_fasta() {
        assert_eq!(OutputFormat::from_extension(Path::new("out.fasta")), OutputFormat::Fasta);
        assert_eq!(OutputFormat::from_extension(Path::new("out.fa")), OutputFormat::Fasta);
        assert_eq!(OutputFormat::from_extension(Path::new("out.fna.gz")), OutputFormat::Fasta);
    }

    #[test]
    fn test_is_gz() {
        assert!(OutputFormat::is_gz(Path::new("file.fastq.gz")));
        assert!(!OutputFormat::is_gz(Path::new("file.fastq")));
    }

    #[test]
    fn test_base_read_id_strips_mate_suffix() {
        assert_eq!(base_read_id(b"HWUSI-EAS100R:6:73:941:1973#0/1"),
                   b"HWUSI-EAS100R:6:73:941:1973#0");
        assert_eq!(base_read_id(b"HWUSI-EAS100R:6:73:941:1973#0/2"),
                   b"HWUSI-EAS100R:6:73:941:1973#0");
        // CASAVA 1.8+ format (no /N suffix, metadata after space)
        assert_eq!(base_read_id(b"HWI-ST1276:71:C1162ACXX:1:1101:1208:2458 1:N:0:CGATGT"),
                   b"HWI-ST1276:71:C1162ACXX:1:1101:1208:2458");
        // SRA split-files variant
        assert_eq!(base_read_id(b"SRR001666.1/1 HWI:1:1:82:1089 length=36"),
                   b"SRR001666.1");
        // No suffix, no description
        assert_eq!(base_read_id(b"read42"), b"read42");
    }
}
