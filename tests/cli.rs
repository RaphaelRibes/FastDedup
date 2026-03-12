use assert_cmd::cargo;
use assert_cmd::Command;
use assert_fs::prelude::*;
use predicates::prelude::*;
use std::fs;

// --------------------------------------------------------
// UTILITAIRE : BUILDER DE COMMANDE
// --------------------------------------------------------

/// Retourne une commande prête à l'emploi pointant vers le binaire fdedup.
/// Centralise le nom du binaire pour faciliter un éventuel renommage.
fn cmd() -> Command {
    cargo::cargo_bin_cmd!("fdedup")
}

// --------------------------------------------------------
// UTILITAIRE : GÉNÉRATION DE SÉQUENCES
// --------------------------------------------------------

const NUCLEOTIDES: [u8; 5] = [b'A', b'C', b'G', b'T', b'N'];

/// Génère une séquence nucléotidique déterministe de 150 bases.
///
/// Effectue exactement 150 étapes de décodage en base 5 sur la valeur
/// `index + offset`. Chaque valeur distincte de `index + offset` produit
/// une séquence distincte (représentation en base 5 sur 150 chiffres,
/// complétée par des 'A' en tête), garantissant l'absence de collision
/// pour tous les indices utilisés dans les tests.
///
/// - `offset = 0`    → séquences R1 / Single-End
/// - `offset = 9999` → séquences R2 (distinctes des R1 pour le même index)
fn generer_sequence(index: usize, offset: usize) -> String {
    const LONGUEUR: usize = 150;
    let mut seq = Vec::with_capacity(LONGUEUR);
    let mut curseur = index.wrapping_add(offset);

    for _ in 0..LONGUEUR {
        seq.push(NUCLEOTIDES[curseur % NUCLEOTIDES.len()]);
        curseur /= NUCLEOTIDES.len();
    }

    String::from_utf8(seq).expect("Séquence UTF-8 invalide")
}

fn generer_fastq(
    dir: &assert_fs::TempDir,
    nom_fichier: &str,
    nb_uniques: usize,
    nb_duplications: usize,
) -> assert_fs::fixture::ChildPath {
    let fichier = dir.child(nom_fichier);
    let mut contenu = String::new();

    for i in 0..nb_uniques {
        let seq = generer_sequence(i, 0);
        let qual: String = (0..150).map(|j| if (i + j) % 10 == 0 { ',' } else { 'F' }).collect();
        contenu.push_str(&format!(
            "@A00123:456:HFWV2DSXX:1:1101:1000:{} 1:N:0:ATGC\n{}\n+\n{}\n",
            1000 + i, seq, qual
        ));
    }

    for i in 0..nb_duplications {
        let seq = generer_sequence(i, 0); // même offset → même séquence → duplication
        let qual: String = (0..150).map(|j| if (i + j) % 10 == 0 { ',' } else { 'F' }).collect();
        contenu.push_str(&format!(
            "@A00123:456:HFWV2DSXX:1:1101:5000:{} 1:N:0:ATGC\n{}\n+\n{}\n",
            1000 + i, seq, qual
        ));
    }

    fichier.write_str(&contenu).unwrap();
    fichier
}

fn generer_fasta(
    dir: &assert_fs::TempDir,
    nom_fichier: &str,
    nb_uniques: usize,
    nb_duplications: usize,
) -> assert_fs::fixture::ChildPath {
    let fichier = dir.child(nom_fichier);
    let mut contenu = String::new();

    for i in 0..nb_uniques {
        contenu.push_str(&format!(
            ">A00123:456:HFWV2DSXX:1:1101:1000:{} 1:N:0:ATGC\n{}\n",
            1000 + i, generer_sequence(i, 0)
        ));
    }

    for i in 0..nb_duplications {
        contenu.push_str(&format!(
            ">A00123:456:HFWV2DSXX:1:1101:5000:{} 1:N:0:ATGC\n{}\n",
            1000 + i, generer_sequence(i, 0)
        ));
    }

    fichier.write_str(&contenu).unwrap();
    fichier
}

fn generer_fastq_paire(
    dir: &assert_fs::TempDir,
    nom_r1: &str,
    nom_r2: &str,
    nb_uniques: usize,
    nb_duplications: usize,
) -> (assert_fs::fixture::ChildPath, assert_fs::fixture::ChildPath) {
    let fichier_r1 = dir.child(nom_r1);
    let fichier_r2 = dir.child(nom_r2);
    let mut contenu_r1 = String::new();
    let mut contenu_r2 = String::new();

    for i in 0..nb_uniques {
        let base_id = format!("@A00123:456:HFWV2DSXX:1:1101:1000:{} ", 1000 + i);
        contenu_r1.push_str(&format!("{}1:N:0:ATGC\n{}\n+\n{}\n", base_id, generer_sequence(i, 0), "F".repeat(150)));
        contenu_r2.push_str(&format!("{}2:N:0:ATGC\n{}\n+\n{}\n", base_id, generer_sequence(i, 9999), "F".repeat(150)));
    }

    for i in 0..nb_duplications {
        let base_id = format!("@A00123:456:HFWV2DSXX:1:1101:5000:{} ", 1000 + i);
        contenu_r1.push_str(&format!("{}1:N:0:ATGC\n{}\n+\n{}\n", base_id, generer_sequence(i, 0), "F".repeat(150)));
        contenu_r2.push_str(&format!("{}2:N:0:ATGC\n{}\n+\n{}\n", base_id, generer_sequence(i, 9999), "F".repeat(150)));
    }

    fichier_r1.write_str(&contenu_r1).unwrap();
    fichier_r2.write_str(&contenu_r2).unwrap();
    (fichier_r1, fichier_r2)
}

fn generer_fasta_paire(
    dir: &assert_fs::TempDir,
    nom_r1: &str,
    nom_r2: &str,
    nb_uniques: usize,
    nb_duplications: usize,
) -> (assert_fs::fixture::ChildPath, assert_fs::fixture::ChildPath) {
    let fichier_r1 = dir.child(nom_r1);
    let fichier_r2 = dir.child(nom_r2);
    let mut contenu_r1 = String::new();
    let mut contenu_r2 = String::new();

    for i in 0..nb_uniques {
        let base_id = format!(">A00123:456:HFWV2DSXX:1:1101:1000:{} ", 1000 + i);
        contenu_r1.push_str(&format!("{}1:N:0:ATGC\n{}\n", base_id, generer_sequence(i, 0)));
        contenu_r2.push_str(&format!("{}2:N:0:ATGC\n{}\n", base_id, generer_sequence(i, 9999)));
    }

    for i in 0..nb_duplications {
        let base_id = format!(">A00123:456:HFWV2DSXX:1:1101:5000:{} ", 1000 + i);
        contenu_r1.push_str(&format!("{}1:N:0:ATGC\n{}\n", base_id, generer_sequence(i, 0)));
        contenu_r2.push_str(&format!("{}2:N:0:ATGC\n{}\n", base_id, generer_sequence(i, 9999)));
    }

    fichier_r1.write_str(&contenu_r1).unwrap();
    fichier_r2.write_str(&contenu_r2).unwrap();
    (fichier_r1, fichier_r2)
}

/// Décompresse un fichier .gz et retourne son contenu sous forme de String.
fn decompresser_gz(chemin: &std::path::Path) -> String {
    let octets = fs::read(chemin).unwrap();
    let mut decoder = flate2::read::GzDecoder::new(&octets[..]);
    let mut contenu = String::new();
    std::io::Read::read_to_string(&mut decoder, &mut contenu).unwrap();
    contenu
}

fn compter_sequences_fastq(contenu: &str) -> usize {
    contenu.lines().filter(|l| l.starts_with('@')).count()
}

fn compter_sequences_fasta(contenu: &str) -> usize {
    contenu.lines().filter(|l| l.starts_with('>')).count()
}

// --------------------------------------------------------
// TESTS SINGLE-END - ERREURS D'ENTRÉE
// --------------------------------------------------------

#[test]
fn test_fichier_entree_inexistant() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    cmd()
        .current_dir(temp_dir.path())
        .arg("-1").arg("fichier_fantome.fastq")
        .assert()
        .failure()
        .stderr(predicate::str::contains("Failed to read input file"));
}

#[test]
fn test_fichier_entree_vide() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let fichier_vide = temp_dir.child("vide.fastq");
    fichier_vide.touch().unwrap();
    cmd()
        .arg("-1").arg(fichier_vide.path())
        .arg("-s")
        .assert()
        .failure()
        .stderr(predicate::str::contains("Failed to read the first two bytes. Is the file empty?"));
}

#[test]
fn test_fichier_entree_mauvais_format() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let fichier_texte = temp_dir.child("mauvais_format.txt");
    fichier_texte.write_str("Ceci n'est pas un fichier FASTQ\nAvec plusieurs lignes.").unwrap();
    cmd()
        .arg("-1").arg(fichier_texte.path())
        .arg("-s")
        .assert()
        .failure()
        .stderr(
            predicate::str::contains("Invalid sequence data")
                .or(predicate::str::contains("Failed to read input file"))
        );
}

#[test]
fn test_fichier_sortie_chemin_inexistant() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let fichier_entree = generer_fastq(&temp_dir, "entree.fastq", 10, 0);
    cmd()
        .arg("-1").arg(fichier_entree.path())
        .arg("-o").arg("/chemin/vers/un/dossier/inexistant/sortie.fastq")
        .assert()
        .failure();
}

#[test]
fn test_fichier_sortie_gz_corrompu() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let fichier_entree = generer_fastq(&temp_dir, "entree.fastq", 10, 0);
    let fichier_sortie_corrompu = temp_dir.child("sortie_corrompue.fastq.gz");
    fichier_sortie_corrompu.write_str("Ceci n'est pas un gzip valide.").unwrap();
    cmd()
        .arg("-1").arg(fichier_entree.path())
        .arg("-o").arg(fichier_sortie_corrompu.path())
        .assert()
        .failure()
        .stderr(predicate::str::contains("Error opening preload file"));
}

// --------------------------------------------------------
// TESTS SINGLE-END - TAUX DE DUPLICATION
// --------------------------------------------------------

#[test]
fn test_taux_duplication_zero_pourcent() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let entree = generer_fastq(&temp_dir, "100_uniques.fastq", 100, 0);
    cmd()
        .arg("-1").arg(entree.path())
        .arg("-s").arg("-v")
        .assert()
        .success()
        .stdout(predicate::str::contains("Processed fragments: 100"))
        .stdout(predicate::str::contains("Duplicates removed: 0.00%"));
}

#[test]
fn test_taux_duplication_9_09_pourcent() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    // 100 uniques + 10 duplicats = 110 lectures, 10/110 ≈ 9.09%
    let entree = generer_fastq(&temp_dir, "110_seq.fastq", 100, 10);
    cmd()
        .arg("-1").arg(entree.path())
        .arg("-s").arg("-v")
        .assert()
        .success()
        .stdout(predicate::str::contains("Processed fragments: 110"))
        .stdout(predicate::str::contains("Duplicates removed: 9.09%"));
}

#[test]
fn test_taux_duplication_50_pourcent() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    // 10 uniques + 10 duplicats = 20 lectures, 10/20 = 50.00%
    let entree = generer_fastq(&temp_dir, "50_pourcent.fastq", 10, 10);
    cmd()
        .arg("-1").arg(entree.path())
        .arg("-s").arg("-v")
        .assert()
        .success()
        .stdout(predicate::str::contains("Processed fragments: 20"))
        .stdout(predicate::str::contains("Duplicates removed: 50.00%"));
}

#[test]
fn test_taux_duplication_entree_toutes_identiques() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    // 10 lectures toutes identiques : la première est unique, les 9 suivantes sont des doublons
    // → 9/10 = 90.00%
    let fichier = temp_dir.child("tous_identiques.fastq");
    let seq = generer_sequence(42, 0);
    let mut contenu = String::new();
    for i in 0..10 {
        contenu.push_str(&format!(
            "@read_{} 1:N:0:ATGC\n{}\n+\n{}\n",
            i, seq, "F".repeat(150)
        ));
    }
    fichier.write_str(&contenu).unwrap();

    cmd()
        .arg("-1").arg(fichier.path())
        .arg("-s").arg("-v")
        .assert()
        .success()
        .stdout(predicate::str::contains("Processed fragments: 10"))
        .stdout(predicate::str::contains("Duplicates removed: 90.00%"));
}

// --------------------------------------------------------
// TESTS SINGLE-END - FORMAT ET COMPRESSION DE SORTIE
// --------------------------------------------------------

#[test]
fn test_sortie_fastq_non_compresse() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let entree = generer_fastq(&temp_dir, "entree.fastq", 10, 0);
    let sortie = temp_dir.child("sortie.fastq");

    cmd()
        .arg("-1").arg(entree.path())
        .arg("-o").arg(sortie.path())
        .assert()
        .success();

    let contenu = fs::read_to_string(sortie.path()).unwrap();
    assert!(contenu.starts_with('@'), "La sortie FASTQ doit commencer par '@'");
    assert!(contenu.contains('+'), "La sortie FASTQ doit contenir la ligne '+'");
    assert_eq!(compter_sequences_fastq(&contenu), 10);
}

#[test]
fn test_sortie_fastq_gz_compresse() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let entree = generer_fastq(&temp_dir, "entree.fastq", 10, 3);
    let sortie = temp_dir.child("sortie.fastq.gz");

    cmd()
        .arg("-1").arg(entree.path())
        .arg("-o").arg(sortie.path())
        .assert()
        .success();

    let octets = fs::read(sortie.path()).unwrap();
    assert_eq!(octets[0], 0x1f, "Magic byte gzip attendu");
    assert_eq!(octets[1], 0x8b, "Magic byte gzip attendu");

    let contenu = decompresser_gz(sortie.path());
    assert!(contenu.starts_with('@'));
    assert_eq!(compter_sequences_fastq(&contenu), 10, "Doit contenir 10 séquences uniques");
}

#[test]
fn test_sortie_fq_extension() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let entree = generer_fastq(&temp_dir, "entree.fastq", 5, 0);
    let sortie = temp_dir.child("sortie.fq");

    cmd()
        .arg("-1").arg(entree.path())
        .arg("-o").arg(sortie.path())
        .assert()
        .success();

    let contenu = fs::read_to_string(sortie.path()).unwrap();
    assert!(contenu.starts_with('@'), "L'extension .fq doit produire du FASTQ");
    assert_eq!(compter_sequences_fastq(&contenu), 5);
}

#[test]
fn test_sortie_fq_gz_compresse() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let entree = generer_fastq(&temp_dir, "entree.fastq", 5, 0);
    let sortie = temp_dir.child("sortie.fq.gz");

    cmd()
        .arg("-1").arg(entree.path())
        .arg("-o").arg(sortie.path())
        .assert()
        .success();

    let octets = fs::read(sortie.path()).unwrap();
    assert_eq!(octets[0], 0x1f);
    assert_eq!(octets[1], 0x8b);

    let contenu = decompresser_gz(sortie.path());
    assert_eq!(compter_sequences_fastq(&contenu), 5);
}

#[test]
fn test_sortie_fasta_extension_fasta() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let entree = generer_fastq(&temp_dir, "entree.fastq", 10, 2);
    let sortie = temp_dir.child("sortie.fasta");

    cmd()
        .arg("-1").arg(entree.path())
        .arg("-o").arg(sortie.path())
        .assert()
        .success();

    let contenu = fs::read_to_string(sortie.path()).unwrap();
    assert!(contenu.starts_with('>'), "La sortie .fasta doit commencer par '>'");
    assert!(!contenu.contains('+'), "La sortie FASTA ne doit pas contenir de ligne '+'");
    assert_eq!(compter_sequences_fasta(&contenu), 10, "Doit contenir 10 séquences uniques");
}

#[test]
fn test_sortie_fasta_extension_fa() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let entree = generer_fastq(&temp_dir, "entree.fastq", 7, 0);
    let sortie = temp_dir.child("sortie.fa");

    cmd()
        .arg("-1").arg(entree.path())
        .arg("-o").arg(sortie.path())
        .assert()
        .success();

    let contenu = fs::read_to_string(sortie.path()).unwrap();
    assert!(contenu.starts_with('>'), "L'extension .fa doit produire du FASTA");
    assert_eq!(compter_sequences_fasta(&contenu), 7);
}

#[test]
fn test_sortie_fasta_extension_fna() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let entree = generer_fastq(&temp_dir, "entree.fastq", 7, 0);
    let sortie = temp_dir.child("sortie.fna");

    cmd()
        .arg("-1").arg(entree.path())
        .arg("-o").arg(sortie.path())
        .assert()
        .success();

    let contenu = fs::read_to_string(sortie.path()).unwrap();
    assert!(contenu.starts_with('>'), "L'extension .fna doit produire du FASTA");
    assert_eq!(compter_sequences_fasta(&contenu), 7);
}

#[test]
fn test_sortie_fasta_gz_compresse() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let entree = generer_fastq(&temp_dir, "entree.fastq", 8, 2);
    let sortie = temp_dir.child("sortie.fasta.gz");

    cmd()
        .arg("-1").arg(entree.path())
        .arg("-o").arg(sortie.path())
        .assert()
        .success();

    let octets = fs::read(sortie.path()).unwrap();
    assert_eq!(octets[0], 0x1f);
    assert_eq!(octets[1], 0x8b);

    let contenu = decompresser_gz(sortie.path());
    assert!(contenu.starts_with('>'));
    assert_eq!(compter_sequences_fasta(&contenu), 8, "Doit contenir 8 séquences uniques");
}

#[test]
fn test_entree_fasta_sortie_fasta() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let entree = generer_fasta(&temp_dir, "entree.fasta", 10, 3);
    let sortie = temp_dir.child("sortie.fasta");

    cmd()
        .arg("-1").arg(entree.path())
        .arg("-o").arg(sortie.path())
        .assert()
        .success();

    let contenu = fs::read_to_string(sortie.path()).unwrap();
    assert!(contenu.starts_with('>'));
    assert_eq!(compter_sequences_fasta(&contenu), 10);
}

/// La conversion FASTA → FASTQ est impossible car les scores de qualité
/// n'existent pas dans les fichiers FASTA. Le binaire doit refuser
/// explicitement cette opération avec un message d'erreur clair.
#[test]
fn test_entree_fasta_sortie_fastq_non_supporte() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let entree = generer_fasta(&temp_dir, "entree.fasta", 10, 0);
    let sortie = temp_dir.child("sortie.fastq");

    cmd()
        .arg("-1").arg(entree.path())
        .arg("-o").arg(sortie.path())
        .assert()
        .failure()
        .stderr(predicate::str::contains("FASTA → FASTQ conversion not supported"));
}

// --------------------------------------------------------
// SINGLE-END TESTS - --force BEHAVIOR
// --------------------------------------------------------

#[test]
fn test_force_overwrites_existing_file() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    // Premier passage : 10 séquences uniques
    let entree_1 = generer_fastq(&temp_dir, "entree_1.fastq", 10, 0);
    let sortie = temp_dir.child("sortie.fastq");

    cmd()
        .arg("-1").arg(entree_1.path())
        .arg("-o").arg(sortie.path())
        .assert()
        .success();

    assert_eq!(compter_sequences_fastq(&fs::read_to_string(sortie.path()).unwrap()), 10);

    // Second pass with --force and a smaller file: must overwrite, not append
    let entree_2 = generer_fastq(&temp_dir, "entree_2.fastq", 4, 0);
    cmd()
        .arg("-1").arg(entree_2.path())
        .arg("-o").arg(sortie.path())
        .arg("--force")
        .assert()
        .success();

    let contenu_final = fs::read_to_string(sortie.path()).unwrap();
    assert_eq!(
        compter_sequences_fastq(&contenu_final),
        4,
        "--force must overwrite the existing file, not append to it"
    );
}

#[test]
fn test_without_force_resumes_from_existing_output() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    // Simule une reprise : la sortie contient déjà 5 séquences,
    // l'entrée contient les mêmes 5 séquences + 5 nouvelles.
    // Without --force, the 5 existing ones should be considered already processed.
    let entree_initiale = generer_fastq(&temp_dir, "entree_initiale.fastq", 5, 0);
    let sortie = temp_dir.child("sortie_reprise.fastq");

    // Premier passage : produit 5 séquences
    cmd()
        .arg("-1").arg(entree_initiale.path())
        .arg("-o").arg(sortie.path())
        .assert()
        .success();

    // Deuxième passage : l'entrée contient les 5 mêmes + 5 nouvelles
    let entree_complete = generer_fastq(&temp_dir, "entree_complete.fastq", 10, 0);
    cmd()
        .arg("-1").arg(entree_complete.path())
        .arg("-o").arg(sortie.path())
        .assert()
        .success();

    let contenu_final = fs::read_to_string(sortie.path()).unwrap();
    assert_eq!(
        compter_sequences_fastq(&contenu_final),
        10,
        "Après reprise, la sortie doit contenir 10 séquences uniques au total"
    );
}

// --------------------------------------------------------
// TESTS SINGLE-END - MODE DE HACHAGE MANUEL
// --------------------------------------------------------

#[test]
fn test_hachage_64_bits_explicite() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let entree = generer_fastq(&temp_dir, "entree.fastq", 20, 5);
    let sortie = temp_dir.child("sortie_64.fastq");

    cmd()
        .arg("-1").arg(entree.path())
        .arg("-o").arg(sortie.path())
        .arg("-H").arg("64")
        .arg("-v")
        .assert()
        .success()
        .stdout(predicate::str::contains("64-bit Hash"));

    let contenu = fs::read_to_string(sortie.path()).unwrap();
    assert_eq!(compter_sequences_fastq(&contenu), 20);
}

#[test]
fn test_hachage_128_bits_explicite() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let entree = generer_fastq(&temp_dir, "entree.fastq", 20, 5);
    let sortie = temp_dir.child("sortie_128.fastq");

    cmd()
        .arg("-1").arg(entree.path())
        .arg("-o").arg(sortie.path())
        .arg("-H").arg("128")
        .arg("-v")
        .assert()
        .success()
        .stdout(predicate::str::contains("128-bit Hash"));

    let contenu = fs::read_to_string(sortie.path()).unwrap();
    assert_eq!(compter_sequences_fastq(&contenu), 20);
}

#[test]
fn test_hachage_64_et_128_bits_produisent_le_meme_resultat() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let entree = generer_fastq(&temp_dir, "entree.fastq", 15, 5);
    let sortie_64 = temp_dir.child("sortie_64.fastq");
    let sortie_128 = temp_dir.child("sortie_128.fastq");

    cmd()
        .arg("-1").arg(entree.path())
        .arg("-o").arg(sortie_64.path())
        .arg("-H").arg("64")
        .assert()
        .success();

    cmd()
        .arg("-1").arg(entree.path())
        .arg("-o").arg(sortie_128.path())
        .arg("-H").arg("128")
        .assert()
        .success();

    let count_64 = compter_sequences_fastq(&fs::read_to_string(sortie_64.path()).unwrap());
    let count_128 = compter_sequences_fastq(&fs::read_to_string(sortie_128.path()).unwrap());
    assert_eq!(count_64, count_128, "Les deux modes de hachage doivent produire le même nombre de séquences uniques");
}

// --------------------------------------------------------
// TESTS SINGLE-END - INTÉGRITÉ APRÈS COMPRESSION
// --------------------------------------------------------

#[test]
fn test_compression_gz_valide_fastq() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let entree = generer_fastq(&temp_dir, "entree.fastq", 50, 10);
    let sortie = temp_dir.child("sortie.fastq.gz");

    cmd()
        .arg("-1").arg(entree.path())
        .arg("-o").arg(sortie.path())
        .assert()
        .success();

    let contenu = decompresser_gz(sortie.path());
    assert!(contenu.starts_with('@'));
    assert_eq!(compter_sequences_fastq(&contenu), 50);
}

#[test]
fn test_compression_gz_valide_fasta() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let entree = generer_fastq(&temp_dir, "entree.fastq", 50, 10);
    let sortie = temp_dir.child("sortie.fasta.gz");

    cmd()
        .arg("-1").arg(entree.path())
        .arg("-o").arg(sortie.path())
        .assert()
        .success();

    let contenu = decompresser_gz(sortie.path());
    assert!(contenu.starts_with('>'));
    assert_eq!(compter_sequences_fasta(&contenu), 50);
}

#[test]
fn test_deduplication_preservee_apres_compression() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let entree = generer_fastq(&temp_dir, "entree.fastq", 10, 10);
    let sortie_gz = temp_dir.child("sortie.fastq.gz");
    let sortie_non_gz = temp_dir.child("sortie.fastq");

    cmd()
        .arg("-1").arg(entree.path())
        .arg("-o").arg(sortie_gz.path())
        .assert()
        .success();

    cmd()
        .arg("-1").arg(entree.path())
        .arg("-o").arg(sortie_non_gz.path())
        .arg("--force")
        .assert()
        .success();

    let count_gz = compter_sequences_fastq(&decompresser_gz(sortie_gz.path()));
    let count_non_gz = compter_sequences_fastq(&fs::read_to_string(sortie_non_gz.path()).unwrap());

    assert_eq!(count_gz, count_non_gz, "Le nombre de séquences doit être identique entre .gz et non-.gz");
    assert_eq!(count_gz, 10, "Doit contenir exactement 10 séquences uniques");
}

// --------------------------------------------------------
// TESTS PAIRED-END - ERREURS ET VALIDATION
// --------------------------------------------------------

#[test]
fn test_paire_erreur_sortie_r2_manquante() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let (r1, r2) = generer_fastq_paire(&temp_dir, "R1.fastq", "R2.fastq", 10, 0);
    let sortie_r1 = temp_dir.child("sortie_r1.fastq"); // chemin dans le temp_dir

    cmd()
        .arg("-1").arg(r1.path())
        .arg("-2").arg(r2.path())
        .arg("-o").arg(sortie_r1.path()) // utilise le temp_dir, pas le répertoire courant
        .assert()
        .failure()
        .stderr(predicate::str::contains("The --output-r2 (-p) argument is required"));
}

#[test]
fn test_paire_desynchronisee_bloquee() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let (r1, r2) = generer_fastq_paire(&temp_dir, "R1_desync.fastq", "R2_desync.fastq", 10, 0);

    // Insérer une lecture supplémentaire en tête de R2 pour créer une désynchronisation
    let mut contenu_r2_corrompu =
        String::from("@A00123:456:HFWV2DSXX:1:1101:9999:9999 2:N:0:ATGC\nATGC\n+\nFFFF\n");
    contenu_r2_corrompu.push_str(&fs::read_to_string(r2.path()).unwrap());
    r2.write_str(&contenu_r2_corrompu).unwrap();

    cmd()
        .arg("-1").arg(r1.path())
        .arg("-2").arg(r2.path())
        .arg("-p").arg(temp_dir.child("dummy_r2.fastq").path())
        .arg("-s")
        .assert()
        .failure()
        .stderr(predicate::str::contains("Critical desynchronization detected"));
}

#[test]
fn test_paire_formats_incoherents() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let (r1, r2) = generer_fastq_paire(&temp_dir, "R1.fastq", "R2.fastq", 10, 0);

    cmd()
        .arg("-1").arg(r1.path())
        .arg("-2").arg(r2.path())
        .arg("-o").arg(temp_dir.child("sortie_R1.fasta").path())
        .arg("-p").arg(temp_dir.child("sortie_R2.fastq").path())
        .assert()
        .failure()
        .stderr(predicate::str::contains("must have the same format"));
}

// --------------------------------------------------------
// TESTS PAIRED-END - TAUX DE DUPLICATION
// --------------------------------------------------------

#[test]
fn test_paire_taux_duplication_33_pourcent() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    // 100 paires uniques + 50 paires dupliquées = 150 paires, 50/150 ≈ 33.33%
    let (r1, r2) = generer_fastq_paire(&temp_dir, "R1.fastq", "R2.fastq", 100, 50);

    cmd()
        .arg("-1").arg(r1.path())
        .arg("-2").arg(r2.path())
        .arg("-p").arg(temp_dir.child("dummy_r2.fastq").path())
        .arg("-s").arg("-v")
        .assert()
        .success()
        .stdout(predicate::str::contains("Paired-End"))
        .stdout(predicate::str::contains("Processed fragments: 150"))
        .stdout(predicate::str::contains("Duplicates removed: 33.33%"));
}

// --------------------------------------------------------
// TESTS PAIRED-END - FORMAT ET COMPRESSION DE SORTIE
// --------------------------------------------------------

#[test]
fn test_paire_sortie_fastq() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let (r1, r2) = generer_fastq_paire(&temp_dir, "R1.fastq", "R2.fastq", 10, 3);
    let sortie_r1 = temp_dir.child("sortie_R1.fastq");
    let sortie_r2 = temp_dir.child("sortie_R2.fastq");

    cmd()
        .arg("-1").arg(r1.path())
        .arg("-2").arg(r2.path())
        .arg("-o").arg(sortie_r1.path())
        .arg("-p").arg(sortie_r2.path())
        .assert()
        .success();

    let contenu_r1 = fs::read_to_string(sortie_r1.path()).unwrap();
    let contenu_r2 = fs::read_to_string(sortie_r2.path()).unwrap();

    assert!(contenu_r1.starts_with('@'));
    assert!(contenu_r2.starts_with('@'));
    assert_eq!(compter_sequences_fastq(&contenu_r1), 10, "R1 doit contenir 10 paires uniques");
    assert_eq!(compter_sequences_fastq(&contenu_r2), 10, "R2 doit contenir 10 paires uniques");
    assert_eq!(
        compter_sequences_fastq(&contenu_r1),
        compter_sequences_fastq(&contenu_r2),
        "R1 et R2 doivent être synchronisés"
    );
}

#[test]
fn test_paire_sortie_fasta() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let (r1, r2) = generer_fastq_paire(&temp_dir, "R1.fastq", "R2.fastq", 10, 3);
    let sortie_r1 = temp_dir.child("sortie_R1.fasta");
    let sortie_r2 = temp_dir.child("sortie_R2.fasta");

    cmd()
        .arg("-1").arg(r1.path())
        .arg("-2").arg(r2.path())
        .arg("-o").arg(sortie_r1.path())
        .arg("-p").arg(sortie_r2.path())
        .assert()
        .success();

    let contenu_r1 = fs::read_to_string(sortie_r1.path()).unwrap();
    let contenu_r2 = fs::read_to_string(sortie_r2.path()).unwrap();

    assert!(contenu_r1.starts_with('>'));
    assert!(contenu_r2.starts_with('>'));
    assert_eq!(compter_sequences_fasta(&contenu_r1), 10);
    assert_eq!(compter_sequences_fasta(&contenu_r2), 10);
}

#[test]
fn test_paire_sortie_fastq_gz() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let (r1, r2) = generer_fastq_paire(&temp_dir, "R1.fastq", "R2.fastq", 10, 3);
    let sortie_r1 = temp_dir.child("sortie_R1.fastq.gz");
    let sortie_r2 = temp_dir.child("sortie_R2.fastq.gz");

    cmd()
        .arg("-1").arg(r1.path())
        .arg("-2").arg(r2.path())
        .arg("-o").arg(sortie_r1.path())
        .arg("-p").arg(sortie_r2.path())
        .assert()
        .success();

    for sortie in [&sortie_r1, &sortie_r2] {
        let octets = fs::read(sortie.path()).unwrap();
        assert_eq!(octets[0], 0x1f);
        assert_eq!(octets[1], 0x8b);
    }

    let contenu_r1 = decompresser_gz(sortie_r1.path());
    let contenu_r2 = decompresser_gz(sortie_r2.path());

    assert!(contenu_r1.starts_with('@'));
    assert_eq!(compter_sequences_fastq(&contenu_r1), 10);
    assert_eq!(compter_sequences_fastq(&contenu_r1), compter_sequences_fastq(&contenu_r2));
}

#[test]
fn test_paire_sortie_fasta_gz() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let (r1, r2) = generer_fastq_paire(&temp_dir, "R1.fastq", "R2.fastq", 10, 3);
    let sortie_r1 = temp_dir.child("sortie_R1.fasta.gz");
    let sortie_r2 = temp_dir.child("sortie_R2.fasta.gz");

    cmd()
        .arg("-1").arg(r1.path())
        .arg("-2").arg(r2.path())
        .arg("-o").arg(sortie_r1.path())
        .arg("-p").arg(sortie_r2.path())
        .assert()
        .success();

    for sortie in [&sortie_r1, &sortie_r2] {
        let octets = fs::read(sortie.path()).unwrap();
        assert_eq!(octets[0], 0x1f);
        assert_eq!(octets[1], 0x8b);
    }

    let contenu_r1 = decompresser_gz(sortie_r1.path());
    let contenu_r2 = decompresser_gz(sortie_r2.path());

    assert!(contenu_r1.starts_with('>'));
    assert_eq!(compter_sequences_fasta(&contenu_r1), 10);
    assert_eq!(compter_sequences_fasta(&contenu_r1), compter_sequences_fasta(&contenu_r2));
}

#[test]
fn test_paire_entree_fasta_sortie_fasta() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let (r1, r2) = generer_fasta_paire(&temp_dir, "R1.fasta", "R2.fasta", 10, 3);
    let sortie_r1 = temp_dir.child("sortie_R1.fasta");
    let sortie_r2 = temp_dir.child("sortie_R2.fasta");

    cmd()
        .arg("-1").arg(r1.path())
        .arg("-2").arg(r2.path())
        .arg("-o").arg(sortie_r1.path())
        .arg("-p").arg(sortie_r2.path())
        .assert()
        .success();

    let contenu_r1 = fs::read_to_string(sortie_r1.path()).unwrap();
    let contenu_r2 = fs::read_to_string(sortie_r2.path()).unwrap();

    assert!(contenu_r1.starts_with('>'));
    assert!(contenu_r2.starts_with('>'));
    assert_eq!(compter_sequences_fasta(&contenu_r1), 10);
    assert_eq!(compter_sequences_fasta(&contenu_r1), compter_sequences_fasta(&contenu_r2));
}