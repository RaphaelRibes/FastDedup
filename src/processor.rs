use std::fs;
use std::fs::{File, OpenOptions};
use std::io::{BufWriter, Write};
use std::path::Path;
use anyhow::{Result, bail, Context};
use flate2::Compression;
use flate2::write::GzEncoder;
use needletail::parse_fastx_file;
use crate::hasher::{HashChecker, SequenceHasher};
use crate::utils::precharger_hashes_existants;

pub(crate) fn executer_deduplication<T: SequenceHasher + 'static>
(
    chemin_entree: &str,
    chemin_sortie: &str,
    force: bool,
    verbose: bool,
    dryrun: bool,
    estimated_capacity: usize
)
    -> Result<(usize, usize)>
{

    let mut checker = HashChecker::<T>::new(estimated_capacity);

    if dryrun {
        if verbose { println!("Option --dryrun activée : calcul du taux de duplication sans écriture de fichier."); }
        let mut reader = parse_fastx_file(chemin_entree).context("Impossible de lire le fichier d'entrée")?;
        let mut sequences_processed = 0;
        let mut dups = 0;

        while let Some(record) = reader.next() {
            let seqrec = record.context("Données de séquence invalides")?;

            let hash = T::hash_seq(&seqrec.seq());
            let is_unique = checker.check(hash);

            if !is_unique {
            } else {
                dups += 1;
            }

            sequences_processed += 1;
        }

        return Ok((sequences_processed, dups));
    }

    // --- PRÉCHARGEMENT ---
    let path_out = Path::new(chemin_sortie);
    let is_gz = path_out.extension().and_then(|s| s.to_str()) == Some("gz");

    if force {
        if verbose { println!("Option --force activée : écrasement de la sortie."); }
    } else {
        let (preloaded, valid_bytes) = precharger_hashes_existants(chemin_sortie, &mut checker, verbose)?;
        if preloaded > 0 && verbose { println!("{} séquences préchargées.", preloaded); }

        if path_out.exists() {
            let current_size = fs::metadata(path_out)?.len();

            if valid_bytes < current_size {
                if is_gz {
                    bail!("Le fichier de sortie (.gz) est corrompu et ne peut pas être tronqué. Utilisez --force pour recommencer.");
                }

                if verbose {
                    println!("Troncature du fichier corrompu de {} à {} octets.", current_size, valid_bytes);
                }
                let file = OpenOptions::new().write(true).open(path_out)?;
                file.set_len(valid_bytes)?;
            }
        }
    }

    // --- PRÉPARATION DE L'ÉCRITURE ---
    let file_out = if path_out.exists() && !force {
        OpenOptions::new().append(true).open(path_out)?
    } else {
        File::create(path_out)?
    };

    let writer: Box<dyn Write> = if is_gz {
        Box::new(GzEncoder::new(file_out, Compression::default()))
    } else {
        Box::new(file_out)
    };
    let mut output = BufWriter::with_capacity(128 * 1024, writer);

    // --- LECTURE ET PARSING ---
    let mut reader = parse_fastx_file(chemin_entree).context("Impossible de lire le fichier d'entrée")?;
    let mut sequences_processed = 0;
    let mut dups = 0;

    while let Some(record) = reader.next() {
        let seqrec = record.context("Données de séquence invalides")?;

        let hash = T::hash_seq(&seqrec.seq());
        let is_unique = checker.check(hash);

        if !is_unique {
            seqrec.write(&mut output, None).context("Erreur lors de l'écriture de la séquence")?;
        } else {
            dups += 1;
        }

        sequences_processed += 1;
    }

    Ok((sequences_processed, dups))
}