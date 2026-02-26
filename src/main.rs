mod cli;
mod hasher;
mod utils;
mod processor;

use std::time::Instant;
use anyhow::{Result, Context};
use crate::hasher::HashType;
use crate::processor::executer_deduplication;
use crate::utils::estimer_capacite_sequences;
use crate::utils::récupérer_la_méthode_de_hachage;
use cli::{Cli, HashMode};
use clap::Parser;

fn dispatch
(
    chemin_entree: &str,
    chemin_sortie: &str,
    estimated_capacity: usize,
    force: bool,
    verbose: bool,
    dryrun: bool,
    hash_type: HashType
) 
    -> Result<(usize, usize)> 
{
    match hash_type {
        HashType::XXH3_64 => {
            if verbose { println!("Mode: Hash 64 bits (Haute performance)"); }
            executer_deduplication::<u64>(chemin_entree, chemin_sortie, force, verbose, dryrun, estimated_capacity)
        },
        HashType::XXH3_128 => {
            if verbose { println!("Mode: Hash 128 bits (Haute sécurité collision)"); }
            executer_deduplication::<u128>(chemin_entree, chemin_sortie, force, verbose, dryrun, estimated_capacity)
        },
    }
}

fn main() -> Result<()> {
    let args = Cli::parse();

    if args.hash.is_some() && args.threshold != 0.01 {
        eprintln!("Warning: --hash specified, the automatic selection threshold ({}) is ignored.", args.threshold);
    }

    if args.verbose {
        println!("Input file: {}", args.input);
        println!("Output file: {}", args.output);
    }

    let cap_input = estimer_capacite_sequences(&args.input)
        .context("Input file not found or inaccessible")?;

    let cap_output = if args.force { 0 } else { estimer_capacite_sequences(&args.output).unwrap_or(0) };
    let total_capacity = cap_input + cap_output;

    let selected_hash_type = match args.hash {
        Some(HashMode::Bit64) => HashType::XXH3_64,
        Some(HashMode::Bit128) => HashType::XXH3_128,
        None => récupérer_la_méthode_de_hachage(total_capacity, args.threshold),
    };

    if args.verbose { println!("Total estimated capacity: {} sequences", total_capacity); }

    let start = Instant::now();

    let (inc, dup) = dispatch(
        &args.input,
        &args.output,
        total_capacity,
        args.force,
        args.verbose,
        args.dryrun,
        selected_hash_type
    )?;

    if args.verbose {
        println!(
            "Processed: {}\nDuplication: {:.2}%",
            inc,
            if inc > 0 { dup as f64 / inc as f64 * 100.0 } else { 0.0 }
        );
        println!("Total execution time: {:.2?}", start.elapsed());
    }

    Ok(())
}