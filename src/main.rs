mod check_hash;

use needletail::parse_fastx_file;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::env;
use std::fs::{self, File, OpenOptions};
use std::io::{self, BufWriter, Write};
use std::path::Path;
use std::time::Instant;
use xxhash_rust::xxh3::xxh3_64;

use crate::check_hash::HashChecker;

fn estimer_capacite_sequences(chemin: &str) -> io::Result<usize> {
    let path = Path::new(chemin);
    if !path.exists() {
        return Ok(0);
    }

    let metadata = fs::metadata(path)?;
    let file_size_bytes = metadata.len();
    let is_gz = path.extension().and_then(|s| s.to_str()) == Some("gz");

    let estimated_capacity = if is_gz {
        (file_size_bytes / 80) as usize
    } else {
        (file_size_bytes / 350) as usize
    };

    Ok(estimated_capacity)
}

fn precharger_hashes_existants(chemin: &str, checker: &mut HashChecker, verbose: bool) -> io::Result<usize> {
    if !Path::new(chemin).exists() {
        return Ok(0);
    }

    if verbose {
        println!("Préchargement des séquences depuis l'output existant...");
    }

    let mut reader = parse_fastx_file(chemin).map_err(|e| io::Error::other(e.to_string()))?;
    let mut count = 0;

    while let Some(record) = reader.next() {
        let seqrec = record.map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e.to_string()))?;
        checker.check(xxh3_64(&seqrec.seq()));
        count += 1;
    }

    Ok(count)
}

fn lire_chaque_quatrieme_ligne(chemin_entree: &str, chemin_sortie: &str, force: bool, verbose: bool) -> io::Result<(usize, usize)> {
    // --- ESTIMATION ET INITIALISATION ---
    let cap_entree = estimer_capacite_sequences(chemin_entree)?;
    let cap_sortie = if force { 0 } else { estimer_capacite_sequences(chemin_sortie).unwrap_or(0) };
    let estimated_capacity = cap_entree + cap_sortie;
    let mut checker = HashChecker::new(estimated_capacity);

    if verbose { println!("Capacité estimée totale : {} séquences", estimated_capacity); }

    // --- PRÉCHARGEMENT ---
    if force {
        if verbose { println!("Option --force activée : écrasement de la sortie."); }
    } else {
        let preloaded = precharger_hashes_existants(chemin_sortie, &mut checker, verbose)?;
        if preloaded > 0 && verbose { println!("{} séquences préchargées.", preloaded); }
    }

    // --- PRÉPARATION DE L'ÉCRITURE ---
    let path_out = Path::new(chemin_sortie);
    let file_out = if path_out.exists() && !force {
        OpenOptions::new().append(true).open(path_out)?
    } else {
        File::create(path_out)?
    };

    let writer: Box<dyn Write> = if path_out.extension().and_then(|s| s.to_str()) == Some("gz") {
        Box::new(GzEncoder::new(file_out, Compression::default()))
    } else {
        Box::new(file_out)
    };
    let mut output = BufWriter::with_capacity(128 * 1024, writer);

    // --- LECTURE ET PARSING ---
    let mut reader = parse_fastx_file(chemin_entree).map_err(|e| io::Error::other(e.to_string()))?;
    let mut sequences_processed = 0;
    let mut dups = 0;

    while let Some(record) = reader.next() {
        let seqrec = record.map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e.to_string()))?;
        let is_unique = !checker.check(xxh3_64(&seqrec.seq()));

        if is_unique {
            seqrec.write(&mut output, None).map_err(|e| io::Error::other(e.to_string()))?;
        } else {
            dups += 1;
        }

        sequences_processed += 1;
    }

    Ok((sequences_processed, dups))
}

fn main() {
    let mut args = env::args();
    let executable_name = args.next().unwrap_or_else(|| "programme".to_string());

    let mut force = false;
    let mut verbose = false;
    let mut positional_args = Vec::new();

    for arg in args {
        if arg == "--force" { force = true; }
        else if arg == "--verbose" || arg == "-v" { verbose = true; }
        else { positional_args.push(arg); }
    }

    if positional_args.is_empty() {
        eprintln!("Usage: {} <fichier_entree> [fichier_sortie] [--force] [--verbose|-v]", executable_name);
        std::process::exit(1);
    }

    let fichier_entree = &positional_args[0];
    let fichier_sortie = if positional_args.len() >= 2 { &positional_args[1] } else { "output.fastq.gz" };

    if verbose {
        println!("Fichier d'entrée : {}", fichier_entree);
        println!("Fichier de sortie : {}", fichier_sortie);
    }

    let debut = Instant::now();
    match lire_chaque_quatrieme_ligne(fichier_entree, fichier_sortie, force, verbose) {
        Ok((inc, dup)) => {
            if verbose {
                println!(
                    "Séquences traitées : {}\n% de duplication : {:.2}%",
                    inc,
                    if inc > 0 { dup as f64 / inc as f64 * 100.0 } else { 0.0 }
                );
            }
        },
        Err(e) => eprintln!("Erreur lors du traitement : {}", e),
    }

    if verbose { println!("Temps d'exécution total : {:.2?}", debut.elapsed()); }
}