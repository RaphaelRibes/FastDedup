use anyhow::{Context, Result};
use needletail::parse_fastx_file;
use std::io::Write;
use std::path::Path;
use std::{fs, io};
use std::fs::{File, OpenOptions};
use flate2::write::GzEncoder;
use flate2::Compression;
use crate::hasher::{SequenceHasher, HashType, HashVerifier};

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

        // Check if it's compressed first
        let without_gz = path_str.strip_suffix(".gz").unwrap_or(&path_str);

        if without_gz.ends_with(".fasta") || without_gz.ends_with(".fa") || without_gz.ends_with(".fna") {
            OutputFormat::Fasta
        } else {
            // Default to FASTQ
            OutputFormat::Fastq
        }
    }

    /// Checks if the file path indicates a gzipped file.
    pub fn is_gz(path: &Path) -> bool {
        path.to_string_lossy().to_lowercase().ends_with(".gz")
    }
}

/// Returns the appropriate `HashType` given the dataset size and collision probability threshold.
pub fn get_hash_method(size: usize, threshold: f64) -> HashType {
    if (2f64 * 2.0f64.powi(64) * threshold).sqrt() < size as f64 {
        HashType::XXH3_128
    } else {
        HashType::XXH3_64
    }
}

/// Prepares a writer for the output file based on the format and whether it should be forced (overwritten) or appended.
pub fn prepare_writer(path: &Path, force: bool) -> Result<(Box<dyn Write>, OutputFormat)> {
    let format = OutputFormat::from_extension(path);
    let is_gz = OutputFormat::is_gz(path);

    let file = if path.exists() && !force {
        OpenOptions::new().append(true).open(path)
            .with_context(|| format!("Failed to open file in append mode: {:?}", path))?
    } else {
        File::create(path)
            .with_context(|| format!("Failed to create file: {:?}", path))?
    };

    let writer: Box<dyn Write> = if is_gz {
        Box::new(GzEncoder::new(file, Compression::default()))
    } else {
        Box::new(file)
    };

    Ok((writer, format))
}

/// Estimates the number of sequences in a file to pre-allocate memory.
pub fn estimate_sequence_capacity<P: AsRef<Path>>(path: P) -> Result<usize> {
    let path_ref = path.as_ref();
    if !path_ref.exists() {
        return Ok(0);
    }

    let metadata = fs::metadata(path_ref)?;
    let file_size_bytes = metadata.len();
    let is_gz = path_ref.extension().and_then(|s| s.to_str()) == Some("gz");

    let estimated_capacity = if is_gz {
        (file_size_bytes / 80) as usize
    } else {
        (file_size_bytes / 350) as usize
    };

    Ok(estimated_capacity)
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

/// Preloads existing hashes into memory from existing output files (Paired-End) and computes valid byte sizes to synchronize truncation.
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

struct ByteCounter(u64);

impl Write for ByteCounter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.0 += buf.len() as u64;
        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}