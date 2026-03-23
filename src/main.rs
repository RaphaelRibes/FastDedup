mod cli;
mod hasher;
mod processor;
mod utils;

use crate::hasher::HashType;
use crate::processor::{execute_deduplication, execute_paired_deduplication};
use crate::utils::estimate_sequence_capacity;
use crate::utils::{get_hash_method, birthday_problem_square_approximation};
use anyhow::{Context, Result};
use clap::Parser;
use cli::{Cli, HashMode};
use std::time::Instant;

/// Dispatches single-end deduplication based on the selected hash type.
fn dispatch(
    input_path: &str,
    output_path: &str,
    estimated_capacity: usize,
    force: bool,
    verbose: bool,
    dry_run: bool,
    compression_level: u32,
    hash_type: &HashType,
) -> Result<(usize, usize)> {
    match hash_type {
        HashType::XXH3_64 => {
            if verbose {
                println!("Single-End Mode: 64-bit Hash");
            }
            execute_deduplication::<u64>(
                input_path,
                output_path,
                force,
                verbose,
                dry_run,
                compression_level,
                estimated_capacity,
            )
        }
        HashType::XXH3_128 => {
            if verbose {
                println!("Single-End Mode: 128-bit Hash");
            }
            execute_deduplication::<u128>(
                input_path,
                output_path,
                force,
                verbose,
                dry_run,
                compression_level,
                estimated_capacity,
            )
        }
    }
}

/// Dispatches paired-end deduplication based on the selected hash type.
fn dispatch_paired(
    input_r1: &str,
    input_r2: &str,
    output_r1: &str,
    output_r2: &str,
    estimated_capacity: usize,
    force: bool,
    verbose: bool,
    dry_run: bool,
    compression_level: u32,
    hash_type: &HashType,
) -> Result<(usize, usize)> {
    match hash_type {
        HashType::XXH3_64 => {
            if verbose {
                println!("Paired-End Mode: 64-bit Combined Hash");
            }
            execute_paired_deduplication::<u64>(
                input_r1,
                input_r2,
                output_r1,
                output_r2,
                force,
                verbose,
                dry_run,
                compression_level,
                estimated_capacity,
            )
        }
        HashType::XXH3_128 => {
            if verbose {
                println!("Paired-End Mode: 128-bit Combined Hash");
            }
            execute_paired_deduplication::<u128>(
                input_r1,
                input_r2,
                output_r1,
                output_r2,
                force,
                verbose,
                dry_run,
                compression_level,
                estimated_capacity,
            )
        }
    }
}

fn main() -> Result<()> {
    let args = Cli::parse();

    if args.hash.is_some() && args.threshold != 0.01 {
        eprintln!(
            "Warning: --hash specifies a hash size, so the automatic selection threshold ({}) is ignored.",
            args.threshold
        );
    }

    if args.verbose {
        println!("Primary input file (R1): {}", args.input);
        println!("Primary output file (R1): {}", args.output);
    }

    // Estimate capacity based only on R1.
    // In Paired-End, 1 pair = 1 fragment = 1 hash, so R1 size is sufficient!
    let cap_input = estimate_sequence_capacity(&args.input)
        .context("Input file not found or inaccessible")?;

    let cap_output = if args.force {
        0
    } else {
        estimate_sequence_capacity(&args.output).unwrap_or(0)
    };

    let total_capacity = cap_input + cap_output;

    let selected_hash_type = match args.hash {
        Some(HashMode::Bit64) => HashType::XXH3_64,
        Some(HashMode::Bit128) => HashType::XXH3_128,
        None => get_hash_method(total_capacity, args.threshold),
    };

    if args.verbose {
        println!("Total estimated hash table capacity: {} fragments", total_capacity);
    }

    let start = Instant::now();

    // Routing Logic: Single-End vs Paired-End
    let (processed, duplicates) = if let Some(input_r2) = &args.input_r2 {

        // Ensure the user provided an output file for R2
        let output_r2 = args.output_r2.as_ref().expect(
            "Critical Error: The --output-r2 (-p) argument is required when --input-r2 (-2) is used."
        );

        if args.verbose {
            println!("Secondary input file (R2): {}", input_r2);
            println!("Secondary output file (R2): {}", output_r2);
            println!("--- Starting Paired-End Processing ---");
        }

        dispatch_paired(
            &args.input,
            input_r2,
            &args.output,
            output_r2,
            total_capacity,
            args.force,
            args.verbose,
            args.dry_run,
            args.compression,
            &selected_hash_type,
        )?
    } else {
        if args.verbose {
            println!("--- Starting Single-End Processing ---");
        }

        dispatch(
            &args.input,
            &args.output,
            total_capacity,
            args.force,
            args.verbose,
            args.dry_run,
            args.compression,
            &selected_hash_type,
        )?
    };

    if args.verbose {
        println!(
            "Processed fragments: {}\nDuplicates removed: {:.2}%",
            processed,
            if processed > 0 {
                duplicates as f64 / processed as f64 * 100.0
            } else {
                0.0
            }
        );
        println!(
            "Estimated collisions: {:.2e}%",
            birthday_problem_square_approximation(processed, &selected_hash_type)*100.0
        );
        println!("Total execution time: {:.2?}", start.elapsed());
    }

    Ok(())
}