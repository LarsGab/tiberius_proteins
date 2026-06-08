#!/usr/bin/env bash
#SBATCH --job-name=vt_setup
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=logs/vt_setup_%j.log

# One-time setup for the vertebrates-test pipeline.
# Equivalent to slurm/00_setup.sh but only downloads what this pipeline needs:
#   - NCBI taxonomy nodes.dmp (for ODB taxon filtering)
#   - OrthoDB Vertebrata partition
#
# Safe to skip if both files already exist on disk.

set -euo pipefail
source "${REPO_DIR}/config_vertebrates_test.sh"

mkdir -p "$WORK_DIR"/{odb/{raw,filtered},peptides,labels,miniprot_labels,miniprot_analysis,logs}

# ── NCBI taxonomy ─────────────────────────────────────────────────────────────
NODES="$WORK_DIR/odb/nodes.dmp"
if [[ -f "$NODES" ]]; then
    echo "[vt_setup] nodes.dmp already present, skipping taxonomy download."
else
    echo "[vt_setup] Downloading NCBI taxonomy dump ..."
    wget -q --show-progress -O "$WORK_DIR/odb/taxdump.tar.gz" "$NCBI_TAXONOMY_URL"
    tar -xzf "$WORK_DIR/odb/taxdump.tar.gz" -C "$WORK_DIR/odb/" nodes.dmp
    rm "$WORK_DIR/odb/taxdump.tar.gz"
    echo "[vt_setup] nodes.dmp extracted."
fi

# ── OrthoDB Vertebrata partition ──────────────────────────────────────────────
ODB="Vertebrata"
DEST="$WORK_DIR/odb/raw/${ODB}.fa.gz"
if [[ -f "$DEST" ]]; then
    echo "[vt_setup] $ODB already on disk, skipping download."
else
    echo "[vt_setup] Downloading ODB partition: $ODB ..."
    wget -q --show-progress -O "$DEST" "${ODB_BASE_URL}/${ODB}.fa.gz"
    echo "[vt_setup] Done: $DEST"
fi

echo "[vt_setup] All downloads complete."
