#!/usr/bin/env bash
#SBATCH --job-name=ins_miniprot
#SBATCH --array=0-9
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=logs/ins_06_miniprot_%A_%a.log

# 10 tasks: 10 insects × 1 level (order).
# Mirrors slurm/06_miniprot.sh including the Diamond pre-filter step.
set -euo pipefail
source /etc/profile.d/modules.sh
module load singularity/3.11.3
source "${REPO_DIR}/config_insects.sh"

TASKS=()
for SPECIES in "${INSECTS_SPECIES_LIST[@]}"; do
    for LEVEL in "${INSECTS_LEVELS[@]}"; do
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
TIBERIUS_PEPTIDES="$WORK_DIR/peptides/${SPECIES}.pep.fa"
OUT_DIR="$WORK_DIR/miniprot/${SPECIES}_excl_${LEVEL_LABEL}"

echo "[ins_miniprot] Species=$SPECIES  Level=$LEVEL_LABEL  Clade=$CLADE"

for f in "$GENOME" "$PROTEINS" "$SCORING_MATRIX" "$TIBERIUS_PEPTIDES"; do
    [[ -f "$f" ]] || { echo "ERROR: file not found: $f"; exit 1; }
done

mkdir -p "$OUT_DIR"

run_tool() {
    if [[ -n "${TIBERIUS_SIF:-}" ]]; then
        singularity exec "$TIBERIUS_SIF" "$@"
    else
        "$@"
    fi
}

# ── Step 0: preprocess protein DB ────────────────────────────────────────────
PROTEINS_FOR_MINIPROT="$PROTEINS"
PREPROCESSED_FA="$OUT_DIR/protein_preprocessed.fa"

if [[ ! -f "$PREPROCESSED_FA" ]]; then
    N_PROT=$(grep -c '^>' "$PROTEINS" || echo 0)
    echo "[ins_miniprot] ODB FASTA has $N_PROT sequences"

    if [[ "$N_PROT" -gt 1000000 ]]; then
        echo "[ins_miniprot] > 1,000,000 proteins – running Diamond soft filter ..."

        DIAMOND_DB="$OUT_DIR/prot_db"
        DIAMOND_HITS="$OUT_DIR/diamond_preprocess.tsv"

        diamond makedb \
            --in "$PROTEINS" \
            --db "$DIAMOND_DB" \
            --threads "$MINIPROT_THREADS"

        diamond blastp \
            --query "$TIBERIUS_PEPTIDES" \
            --db "${DIAMOND_DB}.dmnd" \
            --out "$DIAMOND_HITS" \
            --outfmt 6 qseqid sseqid pident length evalue bitscore qlen slen \
            --evalue 1e-5 \
            --max-target-seqs 200 \
            --very-sensitive \
            --threads "$MINIPROT_THREADS"

        pushd "$OUT_DIR" > /dev/null
        python3 "$TIBERIUS_REPO/tiberius/scripts/rank_species_from_diamond.py" \
            "$DIAMOND_HITS" 13 \
            > "$OUT_DIR/species_rank.tsv"
        popd > /dev/null

        echo "[ins_miniprot] Top species selected:"
        cat "$OUT_DIR/top_species.txt"

        awk '
        BEGIN {
            while ((getline < "'"$OUT_DIR/top_species.txt"'") > 0) {
                wanted[$1] = 1
            }
        }
        /^>/ {
            hdr = substr($0, 2)
            split(hdr, a, /[ \t]/)
            id = a[1]
            species = id
            sub(/_.*/, "", species)
            keep = (species in wanted)
        }
        keep { print }
        ' "$PROTEINS" > "$PREPROCESSED_FA"

        N_FILT=$(grep -c '^>' "$PREPROCESSED_FA" || echo 0)
        echo "[ins_miniprot] Filtered DB: $N_PROT -> $N_FILT sequences"
        PROTEINS_FOR_MINIPROT="$PREPROCESSED_FA"
    else
        echo "[ins_miniprot] <= 1,000,000 proteins – using full DB."
        ln -sf "$PROTEINS" "$PREPROCESSED_FA"
        PROTEINS_FOR_MINIPROT="$PREPROCESSED_FA"
    fi
else
    echo "[ins_miniprot] Preprocessed DB exists, reusing."
    PROTEINS_FOR_MINIPROT="$PREPROCESSED_FA"
fi

echo "[ins_miniprot] Proteins for alignment: $PROTEINS_FOR_MINIPROT"

# ── Step 1: miniprot alignment ────────────────────────────────────────────────
echo "[ins_miniprot] Running miniprot ..."
run_tool miniprot \
    -t "$MINIPROT_THREADS" \
    --aln \
    "$GENOME" \
    "$PROTEINS_FOR_MINIPROT" \
    > "$OUT_DIR/miniprot.aln"

echo "[ins_miniprot] miniprot done: $(wc -l < "$OUT_DIR/miniprot.aln") lines"

# ── Step 2: boundary scorer ───────────────────────────────────────────────────
echo "[ins_miniprot] Running miniprot_boundary_scorer ..."
run_tool miniprot_boundary_scorer \
    -s "$SCORING_MATRIX" \
    -o "$OUT_DIR/miniprot_parsed.gff" \
    < "$OUT_DIR/miniprot.aln"

# ── Step 3: miniprothint ──────────────────────────────────────────────────────
echo "[ins_miniprot] Running miniprothint.py ..."
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

echo "[ins_miniprot] Done."
