#!/usr/bin/env python3
"""Scatter plot of best-hit qcov vs pident, coloured by TP / FP label.

Each dot is one predicted transcript that received at least one Diamond hit.
Transcripts with no hit are excluded from the scatter but their count is shown
in the panel title.

Reads:
  <diamond-dir>/<SPECIES>_excl_<LEVEL>.tsv   (14-col Diamond blastp output)
  <labels-dir>/<SPECIES>_labels.tsv           (transcript_id, class_code)

Writes:
  <out-dir>/scatter_qcov_pident.pdf / .png

Usage:
    python scripts/plot_scatter_qcov_pident.py \
        --diamond-dir $WORK_DIR/diamond/results \
        --labels-dir  $WORK_DIR/labels \
        --out-dir     results/figures/
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── Constants shared with other scripts ──────────────────────────────────────
TP_CLASS_CODES: frozenset[str] = frozenset({"=", "c", "j"})

SPECIES_LABEL = {
    "Arabidopsis_thaliana":      "A. thaliana",
    "Drosophila_melanogaster":   "D. melanogaster",
    "Homo_sapiens":              "H. sapiens",
    "Phaeodactylum_tricornutum": "P. tricornutum",
}
LEVEL_ORDER = ["genus", "order", "subphylum"]

COLOR_TP = "#1f77b4"   # blue
COLOR_FP = "#d62728"   # red

ALPHA = 0.25
DOT_SIZE = 4


# ── I/O ───────────────────────────────────────────────────────────────────────
def load_best_hits(diamond_tsv: Path) -> dict[str, tuple[float, float]]:
    """Return {transcript_id: (pident, qcov)} for the best hit per query."""
    best_bitscore: dict[str, float] = {}
    result: dict[str, tuple[float, float]] = {}
    with open(diamond_tsv) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            c = line.split("\t")
            if len(c) < 14:
                continue
            try:
                qseqid   = c[0]
                pident   = float(c[2])
                qstart   = int(c[6])
                qend     = int(c[7])
                bitscore = float(c[11])
                qlen     = int(c[12])
            except (ValueError, IndexError):
                continue
            qcov = (qend - qstart + 1) / qlen * 100.0 if qlen > 0 else 0.0
            if qseqid not in best_bitscore or bitscore > best_bitscore[qseqid]:
                best_bitscore[qseqid] = bitscore
                result[qseqid] = (pident, qcov)
    return result


def load_labels(labels_tsv: Path) -> dict[str, str]:
    with open(labels_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return {row["transcript_id"]: row["class_code"] for row in reader}


def discover_combos(
    diamond_dir: Path, labels_dir: Path
) -> list[tuple[str, str, Path, Path]]:
    """Return list of (species, level, diamond_tsv, labels_tsv)."""
    pattern = re.compile(r"^(.+)_excl_(genus|order|subphylum)\.tsv$")
    combos = []
    for tsv in sorted(diamond_dir.glob("*_excl_*.tsv")):
        m = pattern.match(tsv.name)
        if not m:
            continue
        species, level = m.group(1), m.group(2)
        labels_tsv = labels_dir / f"{species}_labels.tsv"
        if not labels_tsv.exists():
            print(f"  WARNING: labels not found for {species}, skipping")
            continue
        combos.append((species, level, tsv, labels_tsv))
    return combos


# ── Plotting ──────────────────────────────────────────────────────────────────
def _panel_data(
    hits: dict[str, tuple[float, float]],
    labels: dict[str, str],
    all_ids: set[str],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    """Return (tp_pident, tp_qcov, fp_pident, fp_qcov, n_no_hit)."""
    n_no_hit = sum(1 for t in all_ids if t not in hits)
    tp_pi, tp_qc, fp_pi, fp_qc = [], [], [], []
    for tid in all_ids:
        if tid not in hits:
            continue
        pident, qcov = hits[tid]
        is_tp = labels.get(tid, "u") in TP_CLASS_CODES
        if is_tp:
            tp_pi.append(pident)
            tp_qc.append(qcov)
        else:
            fp_pi.append(pident)
            fp_qc.append(qcov)
    return (np.array(tp_pi), np.array(tp_qc),
            np.array(fp_pi), np.array(fp_qc),
            n_no_hit)


def plot_scatter(
    diamond_dir: Path,
    labels_dir: Path,
    out_dir: Path,
) -> None:
    combos = discover_combos(diamond_dir, labels_dir)
    if not combos:
        raise FileNotFoundError(
            f"No *_excl_(genus|order|subphylum).tsv files in {diamond_dir}"
        )

    # Organise by (species, level) into a 2-D grid
    species_list = sorted({s for s, _, __, ___ in combos},
                          key=lambda s: list(SPECIES_LABEL).index(s)
                          if s in SPECIES_LABEL else s)
    n_species = len(species_list)
    n_levels  = len(LEVEL_ORDER)

    fig, axes = plt.subplots(
        n_species, n_levels,
        figsize=(4.5 * n_levels, 3.8 * n_species),
        sharex=True, sharey=True,
    )
    if n_species == 1:
        axes = [axes]

    # Pre-load labels (one file per species)
    labels_cache: dict[str, dict[str, str]] = {}
    for species in species_list:
        lp = labels_dir / f"{species}_labels.tsv"
        if lp.exists():
            labels_cache[species] = load_labels(lp)

    # Fill panels
    hit_lookup: dict[tuple[str, str], dict[str, tuple[float, float]]] = {}
    for species, level, dtsv, _ in combos:
        print(f"  Loading {species} / {level} ...", flush=True)
        hit_lookup[(species, level)] = load_best_hits(dtsv)

    for row, species in enumerate(species_list):
        labels = labels_cache.get(species, {})
        all_ids = set(labels)
        for col, level in enumerate(LEVEL_ORDER):
            ax = axes[row][col]
            hits = hit_lookup.get((species, level), {})

            tp_pi, tp_qc, fp_pi, fp_qc, n_no_hit = _panel_data(
                hits, labels, all_ids
            )

            # FP underneath TP so TP is visible when they overlap
            ax.scatter(fp_qc, fp_pi, s=DOT_SIZE, c=COLOR_FP,
                       alpha=ALPHA, linewidths=0, rasterized=True, label="FP")
            ax.scatter(tp_qc, tp_pi, s=DOT_SIZE, c=COLOR_TP,
                       alpha=ALPHA, linewidths=0, rasterized=True, label="TP")

            n_tp = len(tp_pi)
            n_fp = len(fp_pi)
            title = (
                f"{SPECIES_LABEL.get(species, species)}  [{level}]\n"
                f"TP={n_tp:,}  FP={n_fp:,}  no-hit={n_no_hit:,}"
            )
            ax.set_title(title, fontsize=8)
            ax.set_xlim(-2, 102)
            ax.set_ylim(-2, 102)
            ax.grid(True, alpha=0.2, linewidth=0.5)

            if col == 0:
                ax.set_ylabel("% identity", fontsize=9)
            if row == n_species - 1:
                ax.set_xlabel("Query coverage (%)", fontsize=9)

    # Shared legend
    tp_patch = mpatches.Patch(color=COLOR_TP, label="TP  (=, c, j)")
    fp_patch = mpatches.Patch(color=COLOR_FP, label="FP  (all other codes)")
    fig.legend(
        handles=[tp_patch, fp_patch],
        loc="lower center", ncol=2, fontsize=9,
        frameon=True, bbox_to_anchor=(0.5, 0.0),
    )

    fig.suptitle(
        "Best Diamond hit: query coverage vs. % identity\n"
        "(coloured by TP / FP; transcripts with no hit excluded)",
        fontsize=12,
    )
    fig.tight_layout(rect=[0, 0.04, 1, 1])

    out_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        path = out_dir / f"scatter_qcov_pident.{ext}"
        fig.savefig(path, dpi=150, bbox_inches="tight")
        print(f"Saved {path}", flush=True)
    plt.close(fig)


# ── Main ──────────────────────────────────────────────────────────────────────
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--diamond-dir", type=Path, required=True,
                    help="Directory with {SPECIES}_excl_{LEVEL}.tsv files")
    ap.add_argument("--labels-dir",  type=Path, required=True,
                    help="Directory with {SPECIES}_labels.tsv files")
    ap.add_argument("--out-dir",     type=Path, default=Path("results/figures"),
                    help="Output directory for figures")
    args = ap.parse_args(argv)

    print(f"Scanning {args.diamond_dir} for Diamond TSVs ...", flush=True)
    plot_scatter(args.diamond_dir, args.labels_dir, args.out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
