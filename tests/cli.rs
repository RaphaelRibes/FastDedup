use assert_cmd::cargo;
use assert_fs::prelude::*;
use predicates::prelude::*;

fn generer_fastq(
    dir: &assert_fs::TempDir,
    nom_fichier: &str,
    nb_uniques: usize,
    nb_duplications: usize,
) -> assert_fs::fixture::ChildPath {
    let fichier = dir.child(nom_fichier);
    let mut contenu = String::new();

    // Pour éviter une sous-estimation de la taille du filtre de Bloom,
    // on allonge la séquence pour que chaque ficher soit assez lourd (environ 500 octets par seq)
    for i in 0..nb_uniques {
        contenu.push_str(&format!("@SEQ_U_{}\n", i));
        let seq = "ACGT".repeat(100) + &format!("{:050}", i);
        contenu.push_str(&format!("{}\n", seq));
        contenu.push_str("+\n");
        contenu.push_str(&"I".repeat(450));
        contenu.push('\n');
    }

    for i in 0..nb_duplications {
        contenu.push_str(&format!("@SEQ_D_{}\n", i));
        let seq = "ACGT".repeat(100) + &format!("{:050}", i); // Copie exacte des premières
        contenu.push_str(&format!("{}\n", seq));
        contenu.push_str("+\n");
        contenu.push_str(&"I".repeat(450));
        contenu.push('\n');
    }

    fichier.write_str(&contenu).unwrap();
    fichier
}

#[test]
fn test_fichier_entree_inexistant() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let mut cmd = cargo::cargo_bin_cmd!("fdedup");
    cmd.current_dir(temp_dir.path())
        .arg("fichier_fantome.fastq")
        .assert()
        .failure()
        .stderr(predicate::str::contains(
            "Impossible de lire le fichier d'entrée",
        ));
}

#[test]
fn test_fichier_entree_vide() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let fichier_vide = temp_dir.child("vide.fastq");
    fichier_vide.touch().unwrap();

    let mut cmd = cargo::cargo_bin_cmd!("fdedup");
    cmd.arg(fichier_vide.path())
        .arg("-s") // simulation
        .assert()
        .failure()
        .stderr(predicate::str::contains(
            "Failed to read the first two bytes. Is the file empty?",
        ));
}

#[test]
fn test_fichier_entree_mauvais_format() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let fichier_texte = temp_dir.child("mauvais_format.txt");
    fichier_texte
        .write_str(
            "Ceci n'est pas un fichier FASTQ mais un fichier texte normal\nAvec plusieurs lignes.",
        )
        .unwrap();

    let mut cmd = cargo::cargo_bin_cmd!("fdedup");
    cmd.arg(fichier_texte.path())
        .arg("-s")
        .assert()
        .failure()
        .stderr(
            predicate::str::contains("Données de séquence invalides")
                .or(predicate::str::contains("Impossible de lire le fichier")),
        );
}

#[test]
fn test_fichier_sortie_sans_permissions_ou_introuvable() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let fichier_entree = generer_fastq(&temp_dir, "entree.fastq", 10, 0);

    let mut cmd = cargo::cargo_bin_cmd!("fdedup");
    cmd.arg(fichier_entree.path())
        .arg("/chemin/vers/un/dossier/qui/n/existe/pas/sortie.fastq")
        .assert()
        .failure();
}

#[test]
fn test_fichier_sortie_gz_corrompu() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let fichier_entree = generer_fastq(&temp_dir, "entree.fastq", 10, 0);

    let fichier_sortie_corrompu = temp_dir.child("sortie_corrompue.fastq.gz");
    fichier_sortie_corrompu
        .write_str("Ceci n'est pas un fichier compressé gzip valide.")
        .unwrap();

    let mut cmd = cargo::cargo_bin_cmd!("fdedup");
    cmd.arg(fichier_entree.path())
        .arg(fichier_sortie_corrompu.path())
        .assert()
        .failure()
        .stderr(predicate::str::contains(
            "Erreur lors de l'ouverture du fichier de préchargement",
        ));
}

#[test]
fn test_taux_duplication_10_et_50_pourcents() {
    let temp_dir = assert_fs::TempDir::new().unwrap();

    // Test 100 séquences uniques
    let seq_100_uniques = generer_fastq(&temp_dir, "100_uniques.fastq", 100, 0);
    let mut cmd = cargo::cargo_bin_cmd!("fdedup");
    cmd.arg(seq_100_uniques.path())
        .arg("-s")
        .arg("-v")
        .assert()
        .success()
        .stdout(predicate::str::contains("Traitées : 100"))
        .stdout(predicate::str::contains("Duplication : 0.00%"));

    // Test 110 séquences (100 uniques + 10 dupliquées)
    let seq_110 = generer_fastq(&temp_dir, "110_seq.fastq", 100, 10);
    let mut cmd2 = cargo::cargo_bin_cmd!("fdedup");
    cmd2.arg(seq_110.path())
        .arg("-s")
        .arg("-v")
        .assert()
        .success()
        .stdout(predicate::str::contains("Traitées : 110"))
        .stdout(predicate::str::contains("Duplication : 9.09%"));

    // Test 150 séquences (100 uniques + 50 dupliquées)
    let seq_150 = generer_fastq(&temp_dir, "150_seq.fastq", 100, 50);
    let mut cmd3 = cargo::cargo_bin_cmd!("fdedup");
    cmd3.arg(seq_150.path())
        .arg("-s")
        .arg("-v")
        .assert()
        .success()
        .stdout(predicate::str::contains("Traitées : 150"))
        .stdout(predicate::str::contains("Duplication : 33.33%"));
}
