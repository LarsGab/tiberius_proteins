#!/usr/bin/env bash
#SBATCH --job-name=vt_miniprot_gffcmp
#SBATCH --array=0-5
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/vt_06_miniprot_gffcmp_%A_%a.log

# 6 tasks: gffcompare Tiberius transcripts vs miniprothint gene models for
# each vertebrate test species at order-level exclusion.

set -euo pipefail
source "${REPO_DIR}/config_vertebrates_test.sh"

TASKS=()
for SPECIES in "${VERTEBRATES_SPECIES_LIST[@]}"; do
    for LEVEL in "${VERTEBRATES_LEVELS[@]}"; do
        TASKS+=("${SPECIES}:${LEVEL}")
    done
done

TASK="${TASKS[$SLURM_ARRAY_TASK_ID]}"
SPECIES="${TASK%%:*}"
LEVEL="${TASK##*:}"
LEVEL_LABEL="${LEVEL_NAMES[$LEVEL]}"
CLADE="${SPECIES_CLADE[$SPECIES]}"

TIBERIUS_GTF="$BENCH_DIR/paper/$CLADE/$SPECIES/results/predictions/tiberius/tiberius_seqlen.gtf"
MINIPROT_DIR="$WORK_DIR/miniprot/${SPECIES}_excl_${LEVEL_LABEL}"
OUT_DIR="$WORK_DIR/miniprot_labels"

echo "[vt_gffcmp] Species=$SPECIES  Level=$LEVEL_LABEL"

for f in "$TIBERIUS_GTF" \
         "$MINIPROT_DIR/miniprot.gtf" \
         "$MINIPROT_DIR/miniprot_trainingGenes.gff"; do
    [[ -f "$f" ]] || { echo "ERROR: file not found: $f"; exit 1; }
done

mkdir -p "$OUT_DIR"

run_gffcmp_and_extract() {
    local ref="$1"
    local label="$2"
    local gc_prefix="$MINIPROT_DIR/gffcmp_${label}"
    local out_tsv="$OUT_DIR/${SPECIES}_excl_${LEVEL_LABEL}_${label}_labels.tsv"

    echo "[vt_gffcmp] Comparing Tiberius vs miniprothint ($label) ..."
    gffcompare \
        -e 3 \
        -T \
        -r "$ref" \
        "$TIBERIUS_GTF" \
        -o "$gc_prefix"

    local annotated_gtf="${gc_prefix}.annotated.gtf"
    [[ -f "$annotated_gtf" ]] || {
        echo "ERROR: gffcompare did not produce $annotated_gtf"
        exit 1
    }

    python3 "$REPO_DIR/scripts/extract_labels.py" \
        "$annotated_gtf" \
        "$out_tsv"

    echo "[vt_gffcmp] Done ($label): $out_tsv"
}

run_gffcmp_and_extract "$MINIPROT_DIR/miniprot.gtf"               "all"
run_gffcmp_and_extract "$MINIPROT_DIR/miniprot_trainingGenes.gff" "training"

echo "[vt_gffcmp] Step complete."
