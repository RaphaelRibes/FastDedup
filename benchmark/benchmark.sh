#!/bin/bash
#SBATCH --job-name=bench_tools
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --partition=cpu-dedicated
#SBATCH --account=dedicated-cpu@cirad-normal
set -e

# Configuration
BIN_DIR="./bin"
FDEDUP="$BIN_DIR/fdedup"
DATA_DIR="./data"
RESULTS_CSV="benchmark_results.csv"
GENOME_FA="$DATA_DIR/hg38.fa"

# Estimation : 150bp SE FASTQ ~ 320 lectures par Mo (environ)
READS_PER_GB=3200000

mkdir -p "$DATA_DIR"

if [ ! -x "$FDEDUP" ]; then
    echo "Erreur : L'exécutable $FDEDUP n'est pas présent ou n'a pas les droits d'exécution."
    exit 1
fi

# 1. Téléchargement et préparation du génome de référence (si nécessaire)
if [ ! -f "$GENOME_FA" ]; then
    echo "Téléchargement du génome humain (GRCh38)..."
    wget -qO- "https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" | gunzip > "$GENOME_FA"
fi

# Initialisation du fichier de résultats
echo "Outil,Taille_Cible_Go,Taille_Reelle_Go,Temp_Ecoule,RAM_Max_Mo" > "$RESULTS_CSV"

# Fonction pour extraire et sauvegarder les métriques
extract_metrics() {
    local tool=$1
    local log_file=$2
    local size_gb=$3
    local actual_size=$4

    local wall_time=$(grep "Elapsed (wall clock) time" "$log_file" | awk '{print $NF}')
    local max_ram_kb=$(grep "Maximum resident set size" "$log_file" | awk '{print $NF}')
    local max_ram_mb=$(echo "scale=2; $max_ram_kb / 1024" | bc)

    echo "${tool},${size_gb},${actual_size},${wall_time},${max_ram_mb}" >> "$RESULTS_CSV"
    echo "-> [$tool] Temps : $wall_time | RAM Max : $max_ram_mb Mo"
}

# 2. Boucle de benchmark
for size_gb in {5,10,15,20,25,30,40,50,60,70,80,90,100}; do

    echo "=================================================="
    echo "Palier ${size_gb}Go : Génération des données PE..."

    # Division par 2 car on génère des paires (R1 + R2 = Taille Cible)
    TARGET_READS=$((size_gb * READS_PER_GB / 2))

    OUT_R1="$DATA_DIR/tmp_${size_gb}GB_R1.fq"
    OUT_R2="$DATA_DIR/tmp_${size_gb}GB_R2.fq"
    OUT_SE="$DATA_DIR/tmp_${size_gb}GB_SE.fq"

    # Génération Paired-End
    wgsim -N $TARGET_READS -1 150 -2 150 -S ${size_gb} "$GENOME_FA" "$OUT_R1" "$OUT_R2" > /dev/null 2>&1

    # Concaténation pour fdedup (pour qu'il traite exactement le même volume de données)
    cat "$OUT_R1" "$OUT_R2" > "$OUT_SE"

    # Mesure de la taille réelle (basée sur le fichier concaténé)
    ACTUAL_SIZE_BYTES=$(wc -c < "$OUT_SE")
    ACTUAL_SIZE_GB=$(echo "scale=2; $ACTUAL_SIZE_BYTES / 1073741824" | bc)
    echo "Taille réelle totale générée : ${ACTUAL_SIZE_GB} Go (Cible: ${size_gb} Go)"

    # ---------------------------------------------------------
    # Benchmarks des outils
    # ---------------------------------------------------------

    # 1. FDedup (sur le fichier concaténé SE)
    echo "Lancement de fdedup..."
    LOG_FDEDUP="$DATA_DIR/time_fdedup.log"
    /usr/bin/env time -v "$FDEDUP" --forcer "$OUT_SE" "$DATA_DIR/out_fdedup.fastq.gz" 2> "$LOG_FDEDUP"
    extract_metrics "fdedup" "$LOG_FDEDUP" "$size_gb" "$ACTUAL_SIZE_GB"

    # 2. Fastp (sur les fichiers PE)
    echo "Lancement de fastp..."
    LOG_FASTP="$DATA_DIR/time_fastp.log"
    /usr/bin/env time -v fastp --in1 "$OUT_R1" --in2 "$OUT_R2" --out1 "$DATA_DIR/out_fastp_R1.fq.gz" --out2 "$DATA_DIR/out_fastp_R2.fq.gz" --dedup --thread 6 2> "$LOG_FASTP" > /dev/null
    extract_metrics "fastp" "$LOG_FASTP" "$size_gb" "$ACTUAL_SIZE_GB"

    # 3. FastUniq (sur les fichiers PE via fichier liste)
    echo "Lancement de fastuniq..."
    LOG_FASTUNIQ="$DATA_DIR/time_fastuniq.log"
    echo -e "$OUT_R1\n$OUT_R2" > "$DATA_DIR/fastuniq_input.txt"
    /usr/bin/env time -v fastuniq -i "$DATA_DIR/fastuniq_input.txt" -t q -o "$DATA_DIR/out_fastuniq_R1.fq" -p "$DATA_DIR/out_fastuniq_R2.fq" 2> "$LOG_FASTUNIQ"
    extract_metrics "fastuniq" "$LOG_FASTUNIQ" "$size_gb" "$ACTUAL_SIZE_GB"

    # 4. Clumpify (sur les fichiers PE, deduplication stricte)
    echo "Lancement de clumpify..."
    LOG_CLUMPIFY="$DATA_DIR/time_clumpify.log"
    /usr/bin/env time -v clumpify.sh in="$OUT_R1" in2="$OUT_R2" out="$DATA_DIR/out_clumpify_R1.fq.gz" out2="$DATA_DIR/out_clumpify_R2.fq.gz" dedupe=t subs=0 2> "$LOG_CLUMPIFY"
    extract_metrics "clumpify" "$LOG_CLUMPIFY" "$size_gb" "$ACTUAL_SIZE_GB"

    # Nettoyage des fichiers lourds de l'itération
    rm -f "$OUT_R1" "$OUT_R2" "$OUT_SE" "$DATA_DIR"/out_* "$DATA_DIR"/time_* "$DATA_DIR"/fastuniq_input.txt
done

echo "=================================================="
echo "Benchmark terminé ! Résultats sauvegardés dans $RESULTS_CSV"