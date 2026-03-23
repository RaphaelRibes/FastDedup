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

/// Prepares a writer for the output file based on the format and whether it should be forced (overwritten) or appended.
pub fn prepare_writer(path: &Path, force: bool, compression_level: u32) -> Result<(Box<dyn Write>, OutputFormat)> {
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
        Box::new(GzEncoder::new(file, Compression::new(compression_level)))
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hasher::HashType;

    // --- birthday_problem_square_approximation ---

    /// À 2^32 séquences avec un hash 64-bit, la probabilité de collision
    /// doit être exactement 0.5 (valeur analytique connue).
    /// Formule : (2^32)^2 / 2^(64+1) = 2^64 / 2^65 = 0.5
    #[test]
    fn test_birthday_64bit_a_2_puissance_32() {
        let n = (u32::MAX as usize) + 1; // 2^32
        let prob = birthday_problem_square_approximation(n, &HashType::XXH3_64);
        assert!(
            (prob - 0.5).abs() < 1e-9,
            "Probabilité attendue 0.5, obtenue {:.12}",
            prob
        );
    }

    /// À 2^32 séquences avec un hash 128-bit, la probabilité doit être
    /// exactement 2^64 / 2^129 = 2^-65 ≈ 2.71e-20.
    #[test]
    fn test_birthday_128bit_a_2_puissance_32() {
        let n = (u32::MAX as usize) + 1; // 2^32
        let prob = birthday_problem_square_approximation(n, &HashType::XXH3_128);
        let expected = (n as f64).powi(2) / 2.0_f64.powi(129);
        assert!(
            (prob - expected).abs() < 1e-30,
            "Probabilité attendue {:.3e}, obtenue {:.3e}",
            expected,
            prob
        );
    }

    /// À 0 séquences, la probabilité de collision doit être 0.
    #[test]
    fn test_birthday_zero_sequences() {
        let prob = birthday_problem_square_approximation(0, &HashType::XXH3_64);
        assert_eq!(prob, 0.0);
    }

    /// À 1 séquence, la probabilité de collision doit être 1 / 2^65 ≈ 2.71e-20.
    #[test]
    fn test_birthday_une_sequence() {
        let prob = birthday_problem_square_approximation(1, &HashType::XXH3_64);
        assert!(prob > 0.0 && prob < 1e-15);
    }

    // --- get_hash_method ---

    /// Une taille très grande (> seuil) avec le threshold par défaut doit
    /// sélectionner XXH3_128.
    #[test]
    fn test_get_hash_method_grand_dataset_selectionne_128() {
        // seuil pour threshold=0.001 : sqrt(2 * 2^64 * 0.001) ≈ 192_059_795
        let grande_taille = 200_000_000_usize;
        let methode = get_hash_method(grande_taille, 0.001);
        assert!(
            matches!(methode, HashType::XXH3_128),
            "Un dataset de {} séquences doit utiliser XXH3_128", grande_taille
        );
    }

    /// Une taille petite (< seuil) avec le threshold par défaut doit
    /// sélectionner XXH3_64.
    #[test]
    fn test_get_hash_method_petit_dataset_selectionne_64() {
        let petite_taille = 1_000_usize;
        let methode = get_hash_method(petite_taille, 0.001);
        assert!(
            matches!(methode, HashType::XXH3_64),
            "Un dataset de {} séquences doit utiliser XXH3_64", petite_taille
        );
    }

    /// Un threshold de 0 (collisions interdites) doit toujours sélectionner
    /// XXH3_128, quelle que soit la taille.
    #[test]
    fn test_get_hash_method_threshold_zero_force_128() {
        // sqrt(0) = 0 < toute taille > 0
        let methode = get_hash_method(1, 0.0);
        assert!(
            matches!(methode, HashType::XXH3_128),
            "Un threshold de 0 doit toujours sélectionner XXH3_128"
        );
    }

    /// Un threshold de 1.0 (toutes les collisions acceptées) doit toujours
    /// sélectionner XXH3_64 pour des tailles pratiques (< 6 milliards).
    #[test]
    fn test_get_hash_method_threshold_un_force_64() {
        // seuil pour threshold=1.0 : sqrt(2 * 2^64) ≈ 6_074_000_999
        let methode = get_hash_method(1_000_000, 1.0);
        assert!(
            matches!(methode, HashType::XXH3_64),
            "Un threshold de 1.0 doit sélectionner XXH3_64 pour un dataset de taille normale"
        );
    }

    /// Vérifie le point de basculement exact : juste en dessous du seuil → 64-bit,
    /// juste au-dessus → 128-bit.
    #[test]
    fn test_get_hash_method_point_de_basculement() {
        // seuil exact pour threshold=0.5 : sqrt(2 * 2^64 * 0.5) = sqrt(2^64) = 2^32 = 4_294_967_296
        let seuil = 4_294_967_296_usize; // 2^32
        let juste_en_dessous = get_hash_method(seuil - 1, 0.5);
        let juste_au_dessus = get_hash_method(seuil + 1, 0.5);
        assert!(
            matches!(juste_en_dessous, HashType::XXH3_64),
            "Juste en dessous du seuil doit être XXH3_64"
        );
        assert!(
            matches!(juste_au_dessus, HashType::XXH3_128),
            "Juste au-dessus du seuil doit être XXH3_128"
        );
    }
}
