mod cli;
mod hasher;
mod processor;
mod utils;
mod writer;
use crate::hasher::HashType;
use crate::processor::{
    DedupConfig, execute_deduplication, execute_deduplication_with_abundance,
    execute_paired_deduplication, execute_paired_deduplication_with_abundance,
};
use crate::utils::{
    birthday_problem_square_approximation, estimate_sequence_capacity, get_hash_method,
};
use anyhow::{Context, Result, bail};
use clap::Parser;
use cli::Cli;
use std::path::Path;
use std::time::Instant;

/// Dispatches single-end deduplication based on the selected hash type.
fn dispatch(
    input_path: &str,
    output_path: &str,
    hash_type: &HashType,
    cfg: &DedupConfig,
) -> Result<(usize, usize)> {
    match hash_type {
        HashType::XXH3_64 => {
            if cfg.verbose {
                println!("Single-End Mode: 64-bit Hash");
            }
            execute_deduplication::<u64>(input_path, output_path, cfg)
        }
        HashType::XXH3_128 => {
            if cfg.verbose {
                println!("Single-End Mode: 128-bit Hash");
            }
            execute_deduplication::<u128>(input_path, output_path, cfg)
        }
    }
}

/// Dispatches single-end abundance-tracking based on the selected hash type.
fn dispatch_abundance(
    input_path: &str,
    output_path: &str,
    hash_type: &HashType,
    cfg: &DedupConfig,
) -> Result<(usize, usize)> {
    match hash_type {
        HashType::XXH3_64 => {
            if cfg.verbose {
                println!("Single-End Mode (Abundance): 64-bit Hash");
            }
            execute_deduplication_with_abundance::<u64>(input_path, output_path, cfg)
        }
        HashType::XXH3_128 => {
            if cfg.verbose {
                println!("Single-End Mode (Abundance): 128-bit Hash");
            }
            execute_deduplication_with_abundance::<u128>(input_path, output_path, cfg)
        }
    }
}

/// Dispatches paired-end deduplication based on the selected hash type.
fn dispatch_paired(
    input_r1: &str,
    input_r2: &str,
    output_r1: &str,
    output_r2: &str,
    hash_type: &HashType,
    cfg: &DedupConfig,
) -> Result<(usize, usize)> {
    match hash_type {
        HashType::XXH3_64 => {
            if cfg.verbose {
                println!("Paired-End Mode: 64-bit Combined Hash");
            }
            execute_paired_deduplication::<u64>(input_r1, input_r2, output_r1, output_r2, cfg)
        }
        HashType::XXH3_128 => {
            if cfg.verbose {
                println!("Paired-End Mode: 128-bit Combined Hash");
            }
            execute_paired_deduplication::<u128>(input_r1, input_r2, output_r1, output_r2, cfg)
        }
    }
}

/// Dispatches paired-end abundance-tracking based on the selected hash type.
fn dispatch_paired_abundance(
    input_r1: &str,
    input_r2: &str,
    output_r1: &str,
    output_r2: &str,
    hash_type: &HashType,
    cfg: &DedupConfig,
) -> Result<(usize, usize)> {
    match hash_type {
        HashType::XXH3_64 => {
            if cfg.verbose {
                println!("Paired-End Mode (Abundance): 64-bit Combined Hash");
            }
            execute_paired_deduplication_with_abundance::<u64>(
                input_r1, input_r2, output_r1, output_r2, cfg,
            )
        }
        HashType::XXH3_128 => {
            if cfg.verbose {
                println!("Paired-End Mode (Abundance): 128-bit Combined Hash");
            }
            execute_paired_deduplication_with_abundance::<u128>(
                input_r1, input_r2, output_r1, output_r2, cfg,
            )
        }
    }
}

fn main() -> Result<()> {
    let args = Cli::parse();

    // --- INPUT VALIDATION ---

    if !Path::new(&args.input).exists() {
        bail!("Input file not found: {}", args.input);
    }

    if let Some(ref r2) = args.input_r2 {
        if !Path::new(r2).exists() {
            bail!("R2 input file not found: {}", r2);
        }
    }

    if args.hash.is_some() && args.threshold != 0.001 {
        eprintln!(
            "Warning: --hash specifies a hash size, so the automatic selection threshold ({}) is ignored.",
            args.threshold
        );
    }

    if args.verbose {
        println!("Primary input file (R1): {}", args.input);
        println!("Primary output file (R1): {}", args.output);
        println!("Expected read length: {} bp", args.read_length);
    }

    let cap_input = estimate_sequence_capacity(&args.input, args.read_length)
        .context("Failed to estimate input file capacity")?;

    let cap_output = if args.force {
        0
    } else {
        estimate_sequence_capacity(&args.output, args.read_length).unwrap_or(0)
    };

    let total_capacity = cap_input + cap_output;

    let selected_hash_type = match args.hash {
        Some(mode) => HashType::from(mode),
        None => get_hash_method(total_capacity, args.threshold),
    };

    if args.verbose {
        println!(
            "Total estimated hash table capacity: {} fragments",
            total_capacity
        );
    }

    if args.dry_run && args.abundance {
        bail!("--abundance is not compatible with --dry-run");
    }

    let cfg = DedupConfig {
        force: args.force,
        verbose: args.verbose,
        dry_run: args.dry_run,
        compression_level: args.compression,
        estimated_capacity: total_capacity,
        read_length: args.read_length,
        abundance_sep: args.abundance_sep.clone(),
    };

    let start = Instant::now();

    let (processed, duplicates) = if let Some(input_r2) = &args.input_r2 {
        let output_r2 = args
            .output_r2
            .as_ref()
            .context("The --output-r2 (-p) argument is required when --input-r2 (-2) is used.")?;

        if args.verbose {
            println!("Secondary input file (R2): {}", input_r2);
            println!("Secondary output file (R2): {}", output_r2);
            println!("--- Starting Paired-End Processing ---");
        }

        if args.abundance {
            dispatch_paired_abundance(
                &args.input,
                input_r2,
                &args.output,
                output_r2,
                &selected_hash_type,
                &cfg,
            )?
        } else {
            dispatch_paired(
                &args.input,
                input_r2,
                &args.output,
                output_r2,
                &selected_hash_type,
                &cfg,
            )?
        }
    } else {
        if args.verbose {
            println!("--- Starting Single-End Processing ---");
        }

        if args.abundance {
            dispatch_abundance(&args.input, &args.output, &selected_hash_type, &cfg)?
        } else {
            dispatch(&args.input, &args.output, &selected_hash_type, &cfg)?
        }
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
            birthday_problem_square_approximation(processed, &selected_hash_type) * 100.0
        );
        println!("Total execution time: {:.2?}", start.elapsed());
    }

    Ok(())
}
