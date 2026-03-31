use fastdedup::hasher::{SequenceHasher, HashVerifier};
use fastdedup::utils::{preload_existing_hashes, preload_existing_paired_hashes, prepare_writer, OutputFormat, compute_write_buffer_size, validate_format_compatibility, write_fasta_record};
use anyhow::{Context, Result, bail};
use needletail::parse_fastx_file;
use std::fs::{OpenOptions, self};
use std::io::{BufWriter, Write};
use std::path::Path;

// ==========================================
// SINGLE-END MODE
// ==========================================

/// Executes deduplication for single-end reads.
pub(crate) fn execute_deduplication<T: SequenceHasher + 'static>(
    input_path: &str,
    output_path: &str,
    force: bool,
    verbose: bool,
    dry_run: bool,
    compression_level: u32,
    estimated_capacity: usize,
    read_length: usize,
) -> Result<(usize, usize)> {
    let mut verifier = HashVerifier::<T>::new(estimated_capacity);

    if dry_run {
        if verbose {
            println!("--dry-run enabled: calculating duplication rate without writing files.");
        }
        let mut reader = parse_fastx_file(input_path).context("Failed to read input file")?;
        let mut processed_sequences = 0;
        let mut duplicates = 0;

        while let Some(record) = reader.next() {
            let seq_record = record.context("Invalid sequence data")?;
            let seq = seq_record.seq();
            let hash = T::hash_sequence(&seq);
            let is_unique = verifier.verify(hash);

            if !is_unique { duplicates += 1; }
            processed_sequences += 1;
        }

        return Ok((processed_sequences, duplicates));
    }

    // --- PRELOADING ---
    let output_path_bound = Path::new(output_path);
    let is_gz = OutputFormat::is_gz(output_path_bound);

    if force {
        if verbose { println!("--force enabled: overwriting output."); }
    } else {
        let (preloads, valid_bytes) = preload_existing_hashes::<T>(output_path, &mut verifier, verbose)?;
        if preloads > 0 && verbose { println!("{} sequences preloaded.", preloads); }

        if output_path_bound.exists() {
            let current_size = fs::metadata(output_path_bound)?.len();

            if valid_bytes < current_size {
                if is_gz { bail!("The output file (.gz) is corrupted and cannot be truncated. Use --force to restart."); }
                if verbose { println!("Truncating corrupted file from {} to {} bytes.", current_size, valid_bytes); }
                let file = OpenOptions::new().write(true).open(output_path_bound)?;
                file.set_len(valid_bytes)?;
            }
        }
    }

    // --- WRITER PREPARATION ---
    let (writer, output_format) = prepare_writer(output_path_bound, force, compression_level)?;
    let buf_size = compute_write_buffer_size(read_length);
    let mut buffered_writer = BufWriter::with_capacity(buf_size, writer);

    if verbose {
        println!("Write buffer size: {} KB", buf_size / 1024);
    }

    // --- READING AND PARSING ---
    let mut reader = parse_fastx_file(input_path).context("Failed to read input file")?;
    let mut processed_sequences = 0;
    let mut duplicates = 0;
    let mut format_validated = false;

    while let Some(record) = reader.next() {
        let seq_record = record.context("Invalid sequence data")?;

        if !format_validated {
            let has_quality = seq_record.qual().is_some_and(|q| !q.is_empty());
            validate_format_compatibility(has_quality, output_format)?;
            format_validated = true;
        }

        let seq = seq_record.seq();
        let hash = T::hash_sequence(&seq);
        let is_unique = verifier.verify(hash);

        if is_unique {
            match output_format {
                OutputFormat::Fasta => {
                    write_fasta_record(&mut buffered_writer, seq_record.id(), &seq)?;
                }
                OutputFormat::Fastq => {
                    seq_record.write(&mut buffered_writer, None)
                        .context("Error writing FASTQ")?;
                }
            }
        } else {
            duplicates += 1;
        }
        processed_sequences += 1;
    }

    buffered_writer.flush().context("Failed to flush output buffer")?;

    Ok((processed_sequences, duplicates))
}

// ==========================================
// PAIRED-END MODE
// ==========================================

/// Executes deduplication for paired-end reads.
pub(crate) fn execute_paired_deduplication<T: SequenceHasher + 'static>(
    input_r1_path: &str,
    input_r2_path: &str,
    output_r1_path: &str,
    output_r2_path: &str,
    force: bool,
    verbose: bool,
    dry_run: bool,
    compression_level: u32,
    estimated_capacity: usize,
    read_length: usize,
) -> Result<(usize, usize)> {
    let mut verifier = HashVerifier::<T>::new(estimated_capacity);

    let mut reader_r1 = parse_fastx_file(input_r1_path).context("Failed to read R1")?;
    let mut reader_r2 = parse_fastx_file(input_r2_path).context("Failed to read R2")?;

    let mut processed_sequences = 0;
    let mut duplicates = 0;

    if dry_run {
        if verbose { println!("--dry-run enabled: calculating Paired-End duplication rate."); }

        while let (Some(record_r1_res), Some(record_r2_res)) = (reader_r1.next(), reader_r2.next()) {
            let seq_r1 = record_r1_res.context("Invalid sequence in R1")?;
            let seq_r2 = record_r2_res.context("Invalid sequence in R2")?;

            let base_id_r1 = seq_r1.id().split(|&b| b == b' ').next().unwrap_or(seq_r1.id());
            let base_id_r2 = seq_r2.id().split(|&b| b == b' ').next().unwrap_or(seq_r2.id());

            if base_id_r1 != base_id_r2 {
                bail!("Critical desynchronization detected! R1: {}, R2: {}",
                    String::from_utf8_lossy(base_id_r1), String::from_utf8_lossy(base_id_r2));
            }

            let r1_seq = seq_r1.seq();
            let r2_seq = seq_r2.seq();
            let combined_hash = T::hash_pair(&r1_seq, &r2_seq);
            if !verifier.verify(combined_hash) { duplicates += 1; }
            processed_sequences += 1;
        }

        return Ok((processed_sequences, duplicates));
    }

    // --- PRELOADING AND SYNCHRONIZED TRUNCATION ---
    let path_r1_bound = Path::new(output_r1_path);
    let path_r2_bound = Path::new(output_r2_path);

    let is_gz_r1 = OutputFormat::is_gz(path_r1_bound);
    let is_gz_r2 = OutputFormat::is_gz(path_r2_bound);

    if force {
        if verbose { println!("--force enabled: overwriting Paired-End outputs."); }
    } else {
        let (preloads, bytes_r1, bytes_r2) = preload_existing_paired_hashes::<T>(
            output_r1_path, output_r2_path, &mut verifier, verbose
        )?;

        if preloads > 0 && verbose {
            println!("{} pairs preloaded and synchronized.", preloads);
        }

        if path_r1_bound.exists() {
            let current_size_r1 = fs::metadata(path_r1_bound)?.len();
            if bytes_r1 < current_size_r1 {
                if is_gz_r1 { bail!("R1 file (.gz) is desynchronized and cannot be truncated. Use --force."); }
                if verbose { println!("Truncating R1 for resynchronization ({} to {} bytes).", current_size_r1, bytes_r1); }
                let file = OpenOptions::new().write(true).open(path_r1_bound)?;
                file.set_len(bytes_r1)?;
            }
        }

        if path_r2_bound.exists() {
            let current_size_r2 = fs::metadata(path_r2_bound)?.len();
            if bytes_r2 < current_size_r2 {
                if is_gz_r2 { bail!("R2 file (.gz) is desynchronized and cannot be truncated. Use --force."); }
                if verbose { println!("Truncating R2 for resynchronization ({} to {} bytes).", current_size_r2, bytes_r2); }
                let file = OpenOptions::new().write(true).open(path_r2_bound)?;
                file.set_len(bytes_r2)?;
            }
        }
    }

    // --- WRITER PREPARATION ---
    let (writer_r1, output_format_r1) = prepare_writer(path_r1_bound, force, compression_level)?;
    let (writer_r2, output_format_r2) = prepare_writer(path_r2_bound, force, compression_level)?;

    if output_format_r1 != output_format_r2 {
        bail!("Output files R1 and R2 must have the same format (FASTA or FASTQ)");
    }

    let buf_size = compute_write_buffer_size(read_length);

    if verbose {
        println!("Write buffer size: {} KB", buf_size / 1024);
    }

    let mut buf_writer_r1 = BufWriter::with_capacity(buf_size, writer_r1);
    let mut buf_writer_r2 = BufWriter::with_capacity(buf_size, writer_r2);

    // --- SYNCHRONIZED MAIN WRITE LOOP ---
    let mut format_validated = false;

    while let (Some(record_r1_res), Some(record_r2_res)) = (reader_r1.next(), reader_r2.next()) {
        let seq_r1 = record_r1_res.context("Invalid sequence in R1")?;
        let seq_r2 = record_r2_res.context("Invalid sequence in R2")?;

        if !format_validated {
            let has_quality = seq_r1.qual().is_some_and(|q| !q.is_empty());
            validate_format_compatibility(has_quality, output_format_r1)?;
            format_validated = true;
        }

        let base_id_r1 = seq_r1.id().split(|&b| b == b' ').next().unwrap_or(seq_r1.id());
        let base_id_r2 = seq_r2.id().split(|&b| b == b' ').next().unwrap_or(seq_r2.id());

        if base_id_r1 != base_id_r2 {
            bail!("Critical desynchronization detected at pair #{}! R1: {}, R2: {}",
                processed_sequences + 1, String::from_utf8_lossy(base_id_r1), String::from_utf8_lossy(base_id_r2));
        }

        let r1_seq = seq_r1.seq();
        let r2_seq = seq_r2.seq();
        let combined_hash = T::hash_pair(&r1_seq, &r2_seq);
        let is_unique = verifier.verify(combined_hash);

        if is_unique {
            match output_format_r1 {
                OutputFormat::Fasta => {
                    write_fasta_record(&mut buf_writer_r1, seq_r1.id(), &r1_seq)?;
                    write_fasta_record(&mut buf_writer_r2, seq_r2.id(), &r2_seq)?;
                }
                OutputFormat::Fastq => {
                    seq_r1.write(&mut buf_writer_r1, None)
                        .context("Error writing R1 FASTQ")?;
                    seq_r2.write(&mut buf_writer_r2, None)
                        .context("Error writing R2 FASTQ")?;
                }
            }
        } else {
            duplicates += 1;
        }

        processed_sequences += 1;
    }

    if reader_r1.next().is_some() || reader_r2.next().is_some() {
        bail!("Desynchronization detected at the end: one file contains more reads than the other!");
    }

    buf_writer_r1.flush().context("Failed to flush R1 output buffer")?;
    buf_writer_r2.flush().context("Failed to flush R2 output buffer")?;

    Ok((processed_sequences, duplicates))
}
