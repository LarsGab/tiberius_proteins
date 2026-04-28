#!/usr/bin/env python3
"""Evaluate miniprothint-based filtering of Tiberius predictions.

For each Tiberius transcript we have two labels:
  1. Ground truth  – class code vs the reference annotation (step 04)
     TP = {=, c, j};  FP = everything else.
  2. Miniprot support – class code vs the miniprothint gene models (step 07).

We sweep several 'support thresholds' (what counts as "supported by miniprot")
and compute precision / recall / F1, exactly as analyze_filter_rules.py does
for Diamond-based filters.

Output TSV: filter_rule | threshold | precision | recall | f1 | n_kept | n_total

Usage:
    analyze_miniprot_filter.py \\
        --annot-labels   $WORK_DIR/labels/{SPECIES}_labels.tsv \\
        --miniprot-labels-all      $WORK_DIR/miniprot_labels/{SPECIES}_excl_{LEVEL}_all_labels.tsv \\
        --miniprot-labels-training $WORK_DIR/miniprot_labels/{SPECIES}_excl_{LEVEL}_training_labels.tsv \\
        --species  SPECIES \\
        --level    LEVEL \\
        --out-dir  $WORK_DIR/miniprot_analysis
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

TP_CLASS_CODES: frozenset[str] = frozenset({"=", "c", "j"})

# Class codes that indicate the Tiberius transcript is supported by a
# miniprothint model (ordered from most to least permissive).
# 'u' = intergenic / no overlap, 'i' = intronic, 'x'/'s' = opposite strand
_NO_SUPPORT_CODES: frozenset[str] = frozenset({"u", "i", "x", "s"})


def load_labels(tsv: Path) -> dict[str, str]:
    with open(tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return {row["transcript_id"]: row["class_code"] for row in reader}


def pr_metrics(
    all_ids: set[str],
    annot_labels: dict[str, str],
    kept: set[str],
) -> tuple[float, float, float]:
    n_tp_total = sum(1 for t in all_ids if annot_labels.get(t, "u") in TP_CLASS_CODES)
    n_tp_kept  = sum(1 for t in kept if annot_labels.get(t, "u") in TP_CLASS_CODES)
    n_fp_kept  = sum(1 for t in kept if annot_labels.get(t, "u") not in TP_CLASS_CODES)

    precision = n_tp_kept / (n_tp_kept + n_fp_kept) if (n_tp_kept + n_fp_kept) > 0 else 0.0
    recall    = n_tp_kept / n_tp_total if n_tp_total > 0 else 0.0
    f1 = (2 * precision * recall / (precision + recall)
          if (precision + recall) > 0 else 0.0)
    return precision, recall, f1


def sweep(
    all_ids: set[str],
    annot_labels: dict[str, str],
    mp_labels_all: dict[str, str],
    mp_labels_training: dict[str, str],
) -> list[dict]:
    rows: list[dict] = []

    def record(rule: str, threshold, kept: set[str]) -> None:
        p, r, f1 = pr_metrics(all_ids, annot_labels, kept)
        rows.append({
            "filter_rule": rule,
            "threshold":   threshold,
            "precision":   f"{p:.4f}",
            "recall":      f"{r:.4f}",
            "f1":          f"{f1:.4f}",
            "n_kept":      len(kept),
            "n_total":     len(all_ids),
        })

    # Baseline
    record("no_filter", 0, all_ids)

    for mp_labels, prefix in [
        (mp_labels_all,      "all"),
        (mp_labels_training, "training"),
    ]:
        # Keep transcripts with any genomic overlap (excludes u/i/x/s)
        kept_any = {
            t for t in all_ids
            if mp_labels.get(t, "u") not in _NO_SUPPORT_CODES
        }
        record(f"{prefix}_any_overlap", 0, kept_any)

        # Keep transcripts with junction-level support {=, c, j}
        kept_junction = {
            t for t in all_ids
            if mp_labels.get(t, "u") in TP_CLASS_CODES
        }
        record(f"{prefix}_junction_support", 0, kept_junction)

        # Keep transcripts with exact or contained match {=, c}
        kept_compat = {
            t for t in all_ids
            if mp_labels.get(t, "u") in {"=", "c"}
        }
        record(f"{prefix}_compatible", 0, kept_compat)

        # Keep transcripts with exact intron-chain match only
        kept_exact = {t for t in all_ids if mp_labels.get(t, "u") == "="}
        record(f"{prefix}_exact_match", 0, kept_exact)

    return rows


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annot-labels",            required=True, type=Path,
                    help="Labels TSV from step 04 (transcript vs reference annotation)")
    ap.add_argument("--miniprot-labels-all",     required=True, type=Path,
                    help="Labels TSV from step 07 (Tiberius vs miniprot.gtf)")
    ap.add_argument("--miniprot-labels-training",required=True, type=Path,
                    help="Labels TSV from step 07 (Tiberius vs miniprot_trainingGenes.gff)")
    ap.add_argument("--species",  required=True)
    ap.add_argument("--level",    required=True)
    ap.add_argument("--out-dir",  required=True, type=Path)
    args = ap.parse_args(argv)

    args.out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[{args.species}/{args.level}] Loading annotation labels ...", flush=True)
    annot_labels = load_labels(args.annot_labels)
    all_ids = set(annot_labels)
    n_tp = sum(1 for t in all_ids if annot_labels[t] in TP_CLASS_CODES)
    print(f"  -> {len(all_ids):,} transcripts  TP={n_tp:,}  FP={len(all_ids)-n_tp:,}", flush=True)

    print(f"[{args.species}/{args.level}] Loading miniprothint support labels ...", flush=True)
    mp_all      = load_labels(args.miniprot_labels_all)
    mp_training = load_labels(args.miniprot_labels_training)
    print(f"  -> all-models    : {sum(1 for t in all_ids if mp_all.get(t,'u') not in _NO_SUPPORT_CODES):,} transcripts with any overlap", flush=True)
    print(f"  -> training-genes: {sum(1 for t in all_ids if mp_training.get(t,'u') not in _NO_SUPPORT_CODES):,} transcripts with any overlap", flush=True)

    rows = sweep(all_ids, annot_labels, mp_all, mp_training)

    out_tsv = args.out_dir / f"{args.species}_{args.level}_miniprot_pr_curves.tsv"
    fieldnames = ["filter_rule", "threshold", "precision", "recall", "f1",
                  "n_kept", "n_total"]
    with open(out_tsv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} rows to {out_tsv}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
