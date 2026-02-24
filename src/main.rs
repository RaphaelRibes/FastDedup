mod check_hash;

use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;
use xxhash_rust::xxh3::xxh3_64;
use std::time::Instant;
use crate::check_hash::HashChecker;
use std::fs;

fn estimer_capacite_sequences(chemin: &str) -> io::Result<usize> {
    let path = Path::new(chemin);
    let metadata = fs::metadata(path)?;
    let file_size_bytes = metadata.len();

    let is_gz = path.extension().and_then(|s| s.to_str()) == Some("gz");

    // Diviseurs heuristiques : ~80 octets/séquence pour un .gz, ~350 octets/séquence en clair.
    let estimated_capacity = if is_gz {
        (file_size_bytes / 80) as usize
    } else {
        (file_size_bytes / 350) as usize
    };

    Ok(estimated_capacity)
}

fn lire_chaque_quatrieme_ligne(chemin_entree: &str, chemin_sortie: &str) -> io::Result<(usize, usize)> {
    // --- LECTURE ---
    let path_in = Path::new(chemin_entree);
    let file_in = File::open(path_in)?;

    let reader: Box<dyn Read> = if path_in.extension().and_then(|s| s.to_str()) == Some("gz") {
        Box::new(MultiGzDecoder::new(file_in))
    } else {
        Box::new(file_in)
    };
    let mut buffered_reader = BufReader::new(reader);

    // --- ÉCRITURE ---
    let path_out = Path::new(chemin_sortie);
    let file_out = File::create(path_out)?;

    let writer: Box<dyn Write> = if path_out.extension().and_then(|s| s.to_str()) == Some("gz") {
        Box::new(GzEncoder::new(file_out, Compression::default()))
    } else {
        Box::new(file_out)
    };

    let mut output = BufWriter::with_capacity(128 * 1024, writer);

    // --- TRAITEMENT ---
    let estimated_capacity = estimer_capacite_sequences(chemin_entree)?;
    println!("Capacité estimée pour le HashChecker : {} séquences", estimated_capacity);
    let mut checker = HashChecker::new(estimated_capacity);

    let mut lines_processed = 0;
    let mut dups = 0;

    // Création des buffers réutilisables
    let mut l1 = Vec::new();
    let mut l2 = Vec::new();
    let mut l3 = Vec::new();
    let mut l4 = Vec::new();

    loop {
        // On vide les buffers avant chaque nouvelle itération (ne libère pas la capacité mémoire)
        l1.clear();
        l2.clear();

        // Lecture de la ligne 1. Si on lit 0 octet, c'est la fin du fichier.
        if buffered_reader.read_until(b'\n', &mut l1)? == 0 {
            break;
        }

        // Lecture de la ligne 2 (la séquence).
        if buffered_reader.read_until(b'\n', &mut l2)? == 0 {
            return Err(io::Error::new(io::ErrorKind::UnexpectedEof, "FASTQ tronqué"));
        }

        // read_until garde le '\n' (et potentiellement le '\r').
        // Pour que le hash soit identique à ta version avec .lines() (qui supprime les retours à la ligne),
        // on prend une slice qui exclut ces caractères.
        let l2_len = l2.len();
        let l2_content = if l2_len > 0 && l2[l2_len - 1] == b'\n' {
            if l2_len > 1 && l2[l2_len - 2] == b'\r' {
                &l2[..l2_len - 2]
            } else {
                &l2[..l2_len - 1]
            }
        } else {
            &l2[..]
        };

        let is_unique = !checker.check(xxh3_64(l2_content));

        l3.clear();
        l4.clear();
        if buffered_reader.read_until(b'\n', &mut l3)? == 0 ||
            buffered_reader.read_until(b'\n', &mut l4)? == 0 {
            return Err(io::Error::new(io::ErrorKind::UnexpectedEof, "FASTQ tronqué"));
        }

        // Ensuite, on décide si on écrit ou si on compte un duplicata
        if is_unique {
            output.write_all(&l1)?;
            output.write_all(&l2)?;
            output.write_all(&l3)?;
            output.write_all(&l4)?;
        } else {
            dups += 1;
        }

        lines_processed += 1;
    }

    Ok((lines_processed, dups))
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
        eprintln!("Usage: {} <fichier_entree> [fichier_sortie]", args[0]);
        std::process::exit(1);
    }

    let fichier_entree = &args[1];

    let fichier_sortie = if args.len() >= 3 {
        &args[2]
    } else {
        "output.fastq.gz"
    };

    println!("Fichier d'entrée : {}", fichier_entree);
    println!("Fichier de sortie : {}", fichier_sortie);

    let debut_sans = Instant::now();
    match lire_chaque_quatrieme_ligne(fichier_entree, fichier_sortie) {
        Ok((inc, dup)) => println!("Nombre de séquences traitées : {}\n%tage de duplication: {}%", inc, dup as f64 / inc as f64 * 100.0),
        Err(e) => eprintln!("Erreur pendant le traitement : {:?}", e),
    }
    let temps_sans = debut_sans.elapsed();

    println!("Temps d'exécution sans multithreading : {:.2?}", temps_sans);
}