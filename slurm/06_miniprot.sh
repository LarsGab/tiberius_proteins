#!/usr/bin/env bash
#SBATCH --job-name=miniprot
#SBATCH --array=0-11
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/06_miniprot_%A_%a.log

# 12 tasks: 4 species × 3 exclusion levels.
# For each task:
#   1. miniprot   – spliced protein-to-genome alignment  -> miniprot.aln
#   2. miniprot_boundary_scorer                         -> miniprot_parsed.gff
#   3. miniprothint.py                                  -> miniprot.gtf
#                                                          miniprot_trainingGenes.gff

set -euo pipefail
source "${REPO_DIR}/config.sh"

TASKS=()
for entry in "${SPECIES_LIST[@]}"; do
    SPECIES="${entry%%:*}"
    for LEVEL in 1 2 3; do
        TASKS+=("${SPECIES}:${LEVEL}")
    done
done

TASK="${TASKS[$SLURM_ARRAY_TASK_ID]}"
SPECIES="${TASK%%:*}"
LEVEL="${TASK##*:}"
LEVEL_LABEL="${LEVEL_NAMES[$LEVEL]}"
CLADE="${SPECIES_CLADE[$SPECIES]}"

GENOME="$BENCH_DIR/$CLADE/$SPECIES/genome.fa"
PROTEINS="$WORK_DIR/odb/filtered/${SPECIES}_excl_${LEVEL_LABEL}.fa"
OUT_DIR="$WORK_DIR/miniprot/${SPECIES}_excl_${LEVEL_LABEL}"

echo "[miniprot] Species=$SPECIES  Level=$LEVEL_LABEL  Clade=$CLADE"
echo "[miniprot] Genome   : $GENOME"
echo "[miniprot] Proteins : $PROTEINS"
echo "[miniprot] Out dir  : $OUT_DIR"

for f in "$GENOME" "$PROTEINS" "$SCORING_MATRIX"; do
    [[ -f "$f" ]] || { echo "ERROR: file not found: $f"; exit 1; }
done

mkdir -p "$OUT_DIR"

# Helper: run with or without Singularity
run_tool() {
    if [[ -n "${TIBERIUS_SIF:-}" ]]; then
        singularity exec "$TIBERIUS_SIF" "$@"
    else
        "$@"
    fi
}

# ── Step 1: miniprot alignment ────────────────────────────────────────────────
echo "[miniprot] Running miniprot ..."
run_tool miniprot \
    -t "$MINIPROT_THREADS" \
    --aln \
    "$GENOME" \
    "$PROTEINS" \
    > "$OUT_DIR/miniprot.aln"

echo "[miniprot] miniprot done: $(wc -l < "$OUT_DIR/miniprot.aln") lines"

# ── Step 2: boundary scorer ───────────────────────────────────────────────────
echo "[miniprot] Running miniprot_boundary_scorer ..."
run_tool miniprot_boundary_scorer \
    -s "$SCORING_MATRIX" \
    -o "$OUT_DIR/miniprot_parsed.gff" \
    < "$OUT_DIR/miniprot.aln"

echo "[miniprot] boundary scorer done: $(wc -l < "$OUT_DIR/miniprot_parsed.gff") GFF lines"

# ── Step 3: miniprothint ──────────────────────────────────────────────────────
echo "[miniprot] Running miniprothint.py ..."
# miniprothint writes output relative to --workdir; we cd there so relative
# paths in its output stay predictable.
pushd "$OUT_DIR" > /dev/null
run_tool miniprothint.py \
    miniprot_parsed.gff \
    --workdir . \
    --ignoreCoverage \
    --topNperSeed 10 \
    --minScoreFraction 0.5
popd > /dev/null

for f in "$OUT_DIR/miniprot.gtf" "$OUT_DIR/miniprot_trainingGenes.gff"; do
    [[ -f "$f" ]] || { echo "ERROR: expected output not found: $f"; exit 1; }
done

echo "[miniprot] Done."
echo "[miniprot]   miniprot.gtf            : $(grep -c $'\t''transcript'$'\t' "$OUT_DIR/miniprot.gtf" || true) transcripts"
echo "[miniprot]   miniprot_trainingGenes  : $(grep -c $'\t''transcript'$'\t' "$OUT_DIR/miniprot_trainingGenes.gff" || true) transcripts"
