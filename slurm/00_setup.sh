#!/usr/bin/env bash
#SBATCH --job-name=prot_setup
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=logs/00_setup_%j.log

set -euo pipefail
source "${REPO_DIR}/config.sh"

mkdir -p "$WORK_DIR"/{odb/{raw,filtered},peptides,diamond/{db,results},labels,analysis,logs}

# ── NCBI taxonomy ─────────────────────────────────────────────────────────────
echo "[setup] Downloading NCBI taxonomy dump ..."
wget -q --show-progress -O "$WORK_DIR/odb/taxdump.tar.gz" "$NCBI_TAXONOMY_URL"
tar -xzf "$WORK_DIR/odb/taxdump.tar.gz" -C "$WORK_DIR/odb/" nodes.dmp
rm "$WORK_DIR/odb/taxdump.tar.gz"
echo "[setup] nodes.dmp extracted."

# ── OrthoDB partitions (download each unique partition once) ──────────────────
declare -A DOWNLOADED

for entry in "${SPECIES_LIST[@]}"; do
    SPECIES="${entry%%:*}"
    ODB="${ODB_PARTITION[$SPECIES]}"
    DEST="$WORK_DIR/odb/raw/${ODB}.fa.gz"

    if [[ -n "${DOWNLOADED[$ODB]:-}" ]]; then
        echo "[setup] $ODB already downloaded, skipping."
        continue
    fi

    if [[ -f "$DEST" ]]; then
        echo "[setup] $ODB already on disk, skipping download."
    else
        echo "[setup] Downloading ODB partition: $ODB ..."
        wget -q --show-progress -O "$DEST" "${ODB_BASE_URL}/${ODB}.fa.gz"
        echo "[setup] Done: $DEST"
    fi
    DOWNLOADED[$ODB]=1
done

echo "[setup] All downloads complete."
