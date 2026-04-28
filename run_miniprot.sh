#!/usr/bin/env bash
# Submit the miniprot-based protein-evidence pipeline to SLURM.
# Assumes steps 01 (ODB filtered FASTAs) and 04 (annotation labels) are done.
#
# Run from the repo root on the cluster:
#   bash run_miniprot.sh

set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export REPO_DIR

source "$REPO_DIR/config.sh"

mkdir -p "$WORK_DIR/logs"

submit() {
    sbatch --parsable "$@"
}

echo "=== Step 6: miniprot + boundary scorer + miniprothint (array 0-11) ==="
JID6=$(submit \
    --partition=pinky,snowball \
    --output="$WORK_DIR/logs/06_miniprot_%A_%a.log" \
    "$REPO_DIR/slurm/06_miniprot.sh")
echo "  Job ID: $JID6"

echo "=== Step 7: gffcompare Tiberius vs miniprothint (array 0-11) ==="
JID7=$(submit \
    --partition=pinky,snowball \
    --dependency=afterok:"$JID6" \
    --output="$WORK_DIR/logs/07_miniprot_gffcmp_%A_%a.log" \
    "$REPO_DIR/slurm/07_miniprot_gffcmp.sh")
echo "  Job ID: $JID7"

echo "=== Step 8: Analyze miniprothint-based filtering (array 0-11) ==="
JID8=$(submit \
    --partition=pinky,snowball \
    --dependency=afterok:"$JID7" \
    --output="$WORK_DIR/logs/08_analyze_miniprot_%A_%a.log" \
    "$REPO_DIR/slurm/08_analyze_miniprot.sh")
echo "  Job ID: $JID8"

echo ""
echo "All jobs submitted."
echo "Monitor with:  squeue -u \$USER"
echo ""
echo "When done, download results:"
echo "  rsync -av brain:$WORK_DIR/miniprot_analysis/ ./results/miniprot_analysis/"
