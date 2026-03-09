use rustc_hash::FxHashSet;
use std::hash::Hash;
use xxhash_rust::xxh3::{xxh3_64, xxh3_128};

/// A hash verifier that checks for sequence uniqueness using a contiguous memory buffer (HashSet).
pub struct HashVerifier<T>
where
    T: Hash + Eq,
{
    memory: FxHashSet<T>,
}

impl<T> HashVerifier<T>
where
    T: Hash + Eq,
{
    /// Initializes the structure with an estimated capacity to minimize reallocations.
    pub fn new(estimated_capacity: usize) -> Self {
        Self {
            memory: FxHashSet::with_capacity_and_hasher(estimated_capacity, Default::default()),
        }
    }

    /// Checks if the hash is already in memory.
    /// Returns `true` if the hash is newly inserted (unique sequence), `false` if it was already present (duplicate).
    pub fn verify(&mut self, hash: T) -> bool {
        self.memory.insert(hash)
    }
}

/// A trait for hashing single and paired sequences across multiple backends.
pub(crate) trait SequenceHasher: Hash + Eq + Copy + Send + Sync {
    /// Hashes a single sequence.
    fn hash_sequence(seq: &[u8]) -> Self;
    /// Hashes a pair of sequences (e.g., paired-end reads).
    fn hash_pair(seq1: &[u8], seq2: &[u8]) -> Self;
}

impl SequenceHasher for u64 {
    #[inline(always)]
    fn hash_sequence(seq: &[u8]) -> Self {
        xxh3_64(seq)
    }

    #[inline(always)]
    fn hash_pair(seq1: &[u8], seq2: &[u8]) -> Self {
        let h1 = xxh3_64(seq1);
        let h2 = xxh3_64(seq2);
        h1 ^ h2.rotate_left(32)
    }
}

impl SequenceHasher for u128 {
    #[inline(always)]
    fn hash_sequence(seq: &[u8]) -> Self {
        xxh3_128(seq)
    }

    #[inline(always)]
    fn hash_pair(seq1: &[u8], seq2: &[u8]) -> Self {
        let h1 = xxh3_128(seq1);
        let h2 = xxh3_128(seq2);
        h1 ^ h2.rotate_left(64)
    }
}

/// The type of hash size to use depending on estimated sequence count.
pub(crate) enum HashType {
    XXH3_64,
    XXH3_128,
}
