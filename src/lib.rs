pub mod hasher;
pub mod utils;

use anyhow::{Context, Result};
use hasher::{HashVerifier, SequenceHasher};
use needletail::parse_fastx_file;
use needletail::parser::SequenceRecord;
use utils::auto_estimate_capacity;

/// A streaming deduplicator that wraps a needletail reader
/// and processes only unique sequences via zero-copy callbacks.
///
/// Instead of returning borrowed records (which is impossible due to
/// needletail's lending-iterator design), this struct drives iteration
/// internally and invokes a user-supplied closure for each unique record.
/// This avoids all per-record allocations while keeping the dedup logic
/// encapsulated and reusable.
///
/// # Examples
///
/// ```no_run
/// use fastdedup::UniqueFastxStream;
///
/// // Without capacity hint (uses a sensible default):
/// let mut stream = UniqueFastxStream::<u64>::new("reads.fastq");
///
/// // With a capacity hint when you know the dataset size:
/// let mut stream = UniqueFastxStream::<u64>::with_capacity("reads.fastq", 50_000_000);
///
/// stream.for_each_unique(|record| {
///     // process unique record…
///     Ok(())
/// }).unwrap();
///
/// let (processed, duplicates) = stream.stats();
/// ```
pub struct UniqueFastxStream<T: SequenceHasher> {
    path: String,
    verifier: HashVerifier<T>,
    processed: usize,
    duplicates: usize,
}

impl<T: SequenceHasher> UniqueFastxStream<T> {
    /// Creates a new stream from a file path, automatically estimating
    /// hash-set capacity from the file size (assumes 150 bp reads).
    ///
    /// If the file does not exist yet or its size cannot be read, a
    /// conservative default is used and the set grows on demand.
    /// For tighter control, use [`with_capacity`](Self::with_capacity).
    pub fn new(path: &str) -> Self {
        let capacity = auto_estimate_capacity(path);
        Self::with_capacity(path, capacity)
    }

    /// Creates a new stream with an explicit capacity hint.
    ///
    /// `estimated_capacity` is forwarded to the internal `HashVerifier`
    /// to pre-allocate the hash set and minimise reallocations.
    /// The estimate does not need to be exact — overshooting wastes a
    /// bit of memory, undershooting just causes a few extra rehashes.
    pub fn with_capacity(path: &str, estimated_capacity: usize) -> Self {
        Self {
            path: path.to_owned(),
            verifier: HashVerifier::new(estimated_capacity),
            processed: 0,
            duplicates: 0,
        }
    }

    /// Drives the reader to completion, calling `func` for every unique
    /// sequence in the file.  Duplicates are silently skipped; parse
    /// errors are propagated immediately.
    ///
    /// The closure receives a zero-copy `&SequenceRecord` whose lifetime
    /// is limited to the callback invocation — no data is copied unless
    /// the caller chooses to do so.
    ///
    /// # Errors
    ///
    /// Returns the first I/O, parse, or callback error encountered.
    #[inline]
    pub fn for_each_unique<F>(&mut self, mut func: F) -> Result<()>
    where
        F: FnMut(&SequenceRecord<'_>) -> Result<()>,
    {
        let mut reader = parse_fastx_file(&self.path)
            .context("Failed to open FASTX file")?;

        while let Some(record_res) = reader.next() {
            let record = record_res.context("Invalid sequence data")?;
            self.processed += 1;

            let hash = T::hash_sequence(&record.seq());
            if self.verifier.verify(hash) {
                func(&record)?;
            } else {
                self.duplicates += 1;
            }
        }

        Ok(())
    }

    /// Returns the current statistics: (processed_sequences, duplicates).
    #[inline]
    pub fn stats(&self) -> (usize, usize) {
        (self.processed, self.duplicates)
    }
}
