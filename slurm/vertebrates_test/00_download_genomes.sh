#!/usr/bin/env bash
#SBATCH --job-name=vt_download
#SBATCH --array=0-5
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/vt_00_download_%A_%a.log

# 6 tasks: one per vertebrate test species.
# Downloads RefSeq genome + annotation via NCBI `datasets` and stages them
# under $BENCH_DIR/Vertebrata/<species>/ so later steps (and the original
# slurm/06_miniprot.sh layout) can find them.
#
# Produces:
#   $BENCH_DIR/Vertebrata/<species>/genome.fa
#   $BENCH_DIR/Vertebrata/<species>/annot_cds.gff
#
# Requires `datasets` (NCBI Datasets CLI) and `unzip` on PATH.

set -euo pipefail
source "${REPO_DIR}/config_vertebrates_test.sh"

SPECIES="${VERTEBRATES_SPECIES_LIST[$SLURM_ARRAY_TASK_ID]}"
CLADE="${SPECIES_CLADE[$SPECIES]}"
ACC="${VERTEBRATES_ACCESSION[$SPECIES]}"

SPECIES_DIR="$BENCH_DIR/$CLADE/$SPECIES"
GENOME="$SPECIES_DIR/genome.fa"
ANNOT="$SPECIES_DIR/annot_cds.gff"

echo "[vt_download] Species=$SPECIES  Clade=$CLADE  Accession=$ACC"
echo "[vt_download] Target dir: $SPECIES_DIR"

if [[ -s "$GENOME" && -s "$ANNOT" ]]; then
    echo "[vt_download] Already present, skipping."
    exit 0
fi

mkdir -p "$SPECIES_DIR"
cd "$SPECIES_DIR"

TMP_ZIP="ncbi.zip"
TMP_DIR="ncbi"

rm -f "$TMP_ZIP"
rm -rf "$TMP_DIR"

echo "[vt_download] Running: datasets download genome accession $ACC"
datasets download genome accession "$ACC" \
    --include genome,gff3 \
    --filename "$TMP_ZIP"

unzip -o -q "$TMP_ZIP" -d "$TMP_DIR"

# NCBI lays out files at: ncbi_dataset/data/<ACC>/{*_genomic.fna,genomic.gff}
DATA_DIR="$TMP_DIR/ncbi_dataset/data/$ACC"
[[ -d "$DATA_DIR" ]] || { echo "ERROR: expected $DATA_DIR"; exit 1; }

# Concatenate any genomic.fna(s) into genome.fa
cat "$DATA_DIR"/*_genomic.fna > "$GENOME"
cp "$DATA_DIR/genomic.gff" "$ANNOT"

# Cleanup downloaded archive
rm -f "$TMP_ZIP"
rm -rf "$TMP_DIR"

[[ -s "$GENOME" ]] || { echo "ERROR: $GENOME empty"; exit 2; }
[[ -s "$ANNOT"  ]] || { echo "ERROR: $ANNOT empty";  exit 2; }

N_SEQ=$(grep -c '^>' "$GENOME" || true)
echo "[vt_download] Done. genome.fa: $N_SEQ sequences;  annot_cds.gff: $(wc -l < "$ANNOT") lines."
