use std::hash::{BuildHasher, Hash};

use anyhow::{Context, Result, bail};
use crossbeam_channel::{Receiver, Sender, bounded};
use needletail::parse_fastx_file;
use rustc_hash::{FxBuildHasher, FxHashSet};

use crate::hasher::SequenceHasher;
use crate::utils::{base_read_id, compute_write_buffer_size, validate_format_compatibility};
use crate::writer::{OutputFormat, prepare_writer, write_fasta_record};
use std::io::{BufWriter, Write};
use std::path::Path;

use super::DedupConfig;

// ==========================================
// TYPES
// ==========================================

/// A thread-safe hash verifier using hash-space sharding to avoid lock contention.
pub(crate) struct ShardedHashVerifier<T: Hash + Eq> {
    shards: Vec<FxHashSet<T>>,
    num_shards: usize,
}

impl<T: Hash + Eq> ShardedHashVerifier<T> {
    pub fn new(capacity: usize, num_shards: usize) -> Self {
        let shard_capacity = (capacity / num_shards).max(1);
        let shards = (0..num_shards)
            .map(|_| FxHashSet::with_capacity_and_hasher(shard_capacity, Default::default()))
            .collect();
        Self { shards, num_shards }
    }

    /// Returns `true` if the hash is newly inserted (unique), `false` if duplicate.
    pub fn verify(&mut self, hash: T) -> bool {
        let idx = self.shard_index(&hash);
        self.shards[idx].insert(hash)
    }

    fn shard_index(&self, hash: &T) -> usize {
        (FxBuildHasher.hash_one(hash) as usize) % self.num_shards
    }
}

/// Owned single-end record data that can be sent across threads.
#[derive(Clone, Debug)]
pub(crate) struct RawRecord {
    pub id: Vec<u8>,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
}

/// Owned paired-end record data.
#[derive(Clone, Debug)]
pub(crate) struct RawPair {
    pub r1: RawRecord,
    pub r2: RawRecord,
}

/// Hashed single-end record, sent from worker to writer.
pub(crate) struct HashedRecord<T> {
    pub hash: T,
    pub record: RawRecord,
}

/// Hashed pair, sent from worker to writer.
pub(crate) struct HashedPair<T> {
    pub hash: T,
    pub pair: RawPair,
}

const CHUNK_SIZE: usize = 4096;

// ==========================================
// SINGLE-END PARALLEL
// ==========================================

/// Executes single-end deduplication using a multi-threaded pipeline.
pub(crate) fn execute_deduplication_parallel<T: SequenceHasher + Send + 'static>(
    input_path: &str,
    output_path: &str,
    cfg: &DedupConfig,
) -> Result<(usize, usize)> {
    let num_workers = resolve_thread_count(cfg.threads);

    let (chunk_tx, chunk_rx): (Sender<Vec<RawRecord>>, Receiver<Vec<RawRecord>>) =
        bounded(num_workers * 2);
    let (hash_tx, hash_rx): (Sender<HashedRecord<T>>, Receiver<HashedRecord<T>>) =
        bounded(num_workers * 8);

    // --- READER THREAD ---
    let reader_handle = std::thread::spawn({
        let chunk_tx = chunk_tx.clone();
        let input_owned = input_path.to_string();
        move || -> Result<(usize, bool)> {
            let mut reader = parse_fastx_file(&input_owned).context("Failed to read input file")?;
            let mut total = 0;
            let mut has_quality = false;
            let mut chunk = Vec::with_capacity(CHUNK_SIZE);

            while let Some(record) = reader.next() {
                let rec = record.context("Invalid sequence")?;
                let id = rec.id().to_vec();
                let seq = rec.seq().to_vec();
                let qual = rec.qual().map(|q| q.to_vec()).unwrap_or_default();

                if !has_quality {
                    has_quality = rec.qual().is_some_and(|q| !q.is_empty());
                }

                chunk.push(RawRecord { id, seq, qual });
                total += 1;

                if chunk.len() >= CHUNK_SIZE {
                    if chunk_tx.send(chunk).is_err() {
                        return Ok((total, has_quality));
                    }
                    chunk = Vec::with_capacity(CHUNK_SIZE);
                }
            }
            if !chunk.is_empty() {
                let _ = chunk_tx.send(chunk);
            }
            Ok((total, has_quality))
        }
    });

    // --- WORKER THREADS ---
    let mut worker_handles = Vec::with_capacity(num_workers);
    for _ in 0..num_workers {
        let chunk_rx = chunk_rx.clone();
        let hash_tx = hash_tx.clone();
        worker_handles.push(std::thread::spawn(move || {
            while let Ok(chunk) = chunk_rx.recv() {
                for record in chunk {
                    let hash = T::hash_sequence(&record.seq);
                    if hash_tx.send(HashedRecord { hash, record }).is_err() {
                        return;
                    }
                }
            }
        }));
    }

    // Drop original senders so channels close when producers finish
    drop(chunk_tx);
    drop(hash_tx);

    // --- WRITER (this thread) ---
    let (_total_reads, has_quality) = reader_handle.join().unwrap()?;

    let output_path = Path::new(output_path);
    let is_gz = OutputFormat::is_gz(output_path);
    let (writer, output_format) = prepare_writer(output_path, true, cfg.compression_level)?;
    let buf_size = compute_write_buffer_size(cfg.read_length);
    let mut buffered_writer = BufWriter::with_capacity(buf_size, writer);

    validate_format_compatibility(has_quality, output_format)?;

    let mut verifier = ShardedHashVerifier::<T>::new(cfg.estimated_capacity, num_workers);
    let mut duplicates = 0;
    let mut processed = 0;

    if cfg.verbose {
        println!(
            "Multi-threaded mode: {} workers, {} shards",
            num_workers, num_workers
        );
    }

    // Receive and process results from workers
    while let Ok(hr) = hash_rx.recv() {
        processed += 1;
        if verifier.verify(hr.hash) {
            match output_format {
                OutputFormat::Fasta => {
                    write_fasta_record(&mut buffered_writer, &hr.record.id, &hr.record.seq)?;
                }
                OutputFormat::Fastq => {
                    buffered_writer.write_all(b"@")?;
                    buffered_writer.write_all(&hr.record.id)?;
                    buffered_writer.write_all(b"\n")?;
                    buffered_writer.write_all(&hr.record.seq)?;
                    buffered_writer.write_all(b"\n+\n")?;
                    buffered_writer.write_all(&hr.record.qual)?;
                    buffered_writer.write_all(b"\n")?;
                }
            }
        } else {
            duplicates += 1;
        }
    }

    // Verify all output is flushed
    if processed > 0 {
        buffered_writer.flush().context("Failed to flush output")?;
    }

    // Worker handles must finish since hash_rx closed
    for h in worker_handles {
        h.join().unwrap();
    }

    if cfg.verbose && is_gz {
        println!(
            "Output file written (gzip compressed): {}",
            output_path.display()
        );
    }

    Ok((processed, duplicates))
}

// ==========================================
// PAIRED-END PARALLEL
// ==========================================

/// Executes paired-end deduplication using a multi-threaded pipeline.
pub(crate) fn execute_paired_deduplication_parallel<T: SequenceHasher + Send + 'static>(
    input_r1_path: &str,
    input_r2_path: &str,
    output_r1_path: &str,
    output_r2_path: &str,
    cfg: &DedupConfig,
) -> Result<(usize, usize)> {
    let num_workers = resolve_thread_count(cfg.threads);

    let (chunk_tx, chunk_rx): (Sender<Vec<RawPair>>, Receiver<Vec<RawPair>>) =
        bounded(num_workers * 2);
    let (hash_tx, hash_rx): (Sender<HashedPair<T>>, Receiver<HashedPair<T>>) =
        bounded(num_workers * 8);

    // --- READER THREAD ---
    let reader_handle = std::thread::spawn({
        let chunk_tx = chunk_tx.clone();
        let r1_owned = input_r1_path.to_string();
        let r2_owned = input_r2_path.to_string();
        move || -> Result<(usize, bool)> {
            let mut reader_r1 = parse_fastx_file(&r1_owned).context("Failed to read R1")?;
            let mut reader_r2 = parse_fastx_file(&r2_owned).context("Failed to read R2")?;

            let mut total = 0;
            let mut has_quality = false;
            let mut chunk = Vec::with_capacity(CHUNK_SIZE);

            while let Some(r1) = reader_r1.next() {
                let r1 = r1.context("Invalid R1 sequence")?;
                let r2 = match reader_r2.next() {
                    Some(r) => r.context("Invalid R2 sequence")?,
                    None => {
                        bail!("Desynchronization: R1 has more records than R2");
                    }
                };

                // Validate pair synchronization
                let base_r1 = base_read_id(r1.id()).to_vec();
                let base_r2 = base_read_id(r2.id()).to_vec();
                if base_r1 != base_r2 {
                    bail!(
                        "Desynchronization: R1 base ID {:?} != R2 base ID {:?}",
                        String::from_utf8_lossy(&base_r1),
                        String::from_utf8_lossy(&base_r2)
                    );
                }

                if !has_quality {
                    has_quality = r1.qual().is_some_and(|q| !q.is_empty());
                }

                let pair = RawPair {
                    r1: RawRecord {
                        id: r1.id().to_vec(),
                        seq: r1.seq().to_vec(),
                        qual: r1.qual().map(|q| q.to_vec()).unwrap_or_default(),
                    },
                    r2: RawRecord {
                        id: r2.id().to_vec(),
                        seq: r2.seq().to_vec(),
                        qual: r2.qual().map(|q| q.to_vec()).unwrap_or_default(),
                    },
                };

                chunk.push(pair);
                total += 1;

                if chunk.len() >= CHUNK_SIZE {
                    if chunk_tx.send(chunk).is_err() {
                        return Ok((total, has_quality));
                    }
                    chunk = Vec::with_capacity(CHUNK_SIZE);
                }
            }

            // Check for R2 having extra records
            if reader_r2.next().is_some() {
                bail!("Desynchronization: R2 has more records than R1");
            }

            if !chunk.is_empty() {
                let _ = chunk_tx.send(chunk);
            }
            Ok((total, has_quality))
        }
    });

    // --- WORKER THREADS ---
    let mut worker_handles = Vec::with_capacity(num_workers);
    for _ in 0..num_workers {
        let chunk_rx = chunk_rx.clone();
        let hash_tx = hash_tx.clone();
        worker_handles.push(std::thread::spawn(move || {
            while let Ok(chunk) = chunk_rx.recv() {
                for pair in chunk {
                    let hash = T::hash_pair(&pair.r1.seq, &pair.r2.seq);
                    if hash_tx.send(HashedPair { hash, pair }).is_err() {
                        return;
                    }
                }
            }
        }));
    }

    drop(chunk_tx);
    drop(hash_tx);

    // --- WRITER (this thread) ---
    let (_total_reads, has_quality) = match reader_handle.join() {
        Ok(Ok(result)) => result,
        Ok(Err(e)) => return Err(e),
        Err(e) => bail!("Reader thread panicked: {:?}", e),
    };

    let path_r1 = Path::new(output_r1_path);
    let path_r2 = Path::new(output_r2_path);
    let (writer_r1, output_format) = prepare_writer(path_r1, true, cfg.compression_level)?;
    let (writer_r2, output_format_r2) = prepare_writer(path_r2, true, cfg.compression_level)?;

    if output_format != output_format_r2 {
        bail!("Output files R1 and R2 must have the same format");
    }

    let buf_size = compute_write_buffer_size(cfg.read_length);
    let mut buf_r1 = BufWriter::with_capacity(buf_size, writer_r1);
    let mut buf_r2 = BufWriter::with_capacity(buf_size, writer_r2);

    validate_format_compatibility(has_quality, output_format)?;

    let mut verifier = ShardedHashVerifier::<T>::new(cfg.estimated_capacity, num_workers);
    let mut duplicates = 0;
    let mut processed = 0;

    if cfg.verbose {
        println!(
            "Multi-threaded PE mode: {} workers, {} shards",
            num_workers, num_workers
        );
    }

    while let Ok(hp) = hash_rx.recv() {
        processed += 1;
        if verifier.verify(hp.hash) {
            match output_format {
                OutputFormat::Fasta => {
                    write_fasta_record(&mut buf_r1, &hp.pair.r1.id, &hp.pair.r1.seq)?;
                    write_fasta_record(&mut buf_r2, &hp.pair.r2.id, &hp.pair.r2.seq)?;
                }
                OutputFormat::Fastq => {
                    buf_r1.write_all(b"@")?;
                    buf_r1.write_all(&hp.pair.r1.id)?;
                    buf_r1.write_all(b"\n")?;
                    buf_r1.write_all(&hp.pair.r1.seq)?;
                    buf_r1.write_all(b"\n+\n")?;
                    buf_r1.write_all(&hp.pair.r1.qual)?;
                    buf_r1.write_all(b"\n")?;

                    buf_r2.write_all(b"@")?;
                    buf_r2.write_all(&hp.pair.r2.id)?;
                    buf_r2.write_all(b"\n")?;
                    buf_r2.write_all(&hp.pair.r2.seq)?;
                    buf_r2.write_all(b"\n+\n")?;
                    buf_r2.write_all(&hp.pair.r2.qual)?;
                    buf_r2.write_all(b"\n")?;
                }
            }
        } else {
            duplicates += 1;
        }
    }

    if processed > 0 {
        buf_r1.flush()?;
        buf_r2.flush()?;
    }

    for h in worker_handles {
        h.join().unwrap();
    }

    Ok((processed, duplicates))
}

// ==========================================
// HELPERS
// ==========================================

fn resolve_thread_count(requested: usize) -> usize {
    if requested == 0 {
        std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(4)
    } else {
        requested
    }
}
