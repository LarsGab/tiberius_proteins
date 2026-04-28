#!/usr/bin/env python3
"""Sweep Diamond-based filter rules and compute precision / recall tables.

For a given (species, ODB exclusion level) pair this script:
  1. Reads the Diamond blastp TSV (must include qlen as col 13, slen as col 14).
  2. Reads the per-transcript label TSV from extract_labels.py.
  3. Computes the best hit per query transcript (by bitscore).
  4. Sweeps a family of filter rules and records precision / recall / F1.

A prediction is a TRUE POSITIVE when its gffcompare class_code is in
TP_CLASS_CODES.  All other predictions are FALSE POSITIVES.

A filter rule keeps a transcript if it has a Diamond hit that passes the
given criterion.  Transcripts with NO hit are always removed by the filter.
'no_filter' (baseline) keeps every transcript regardless.

Output: TSV  filter_rule | threshold | precision | recall | f1 | n_kept | n_total

Usage:
    analyze_filter_rules.py --diamond TSV --labels TSV
                            --species NAME --level LABEL --out-dir DIR
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import NamedTuple

# Class codes counted as true positives.
# '='  exact intron-chain match (strictest)
# 'c'  contained in reference, intron-compatible
# 'j'  potential novel isoform with ≥1 matching junction
TP_CLASS_CODES: frozenset[str] = frozenset({"=", "c", "j"})


class Hit(NamedTuple):
    pident:    float   # % identity
    query_cov: float   # (qend-qstart+1) / qlen * 100
    evalue:    float
    bitscore:  float


# ── I/O helpers ───────────────────────────────────────────────────────────────

def load_best_hits(diamond_tsv: Path) -> dict[str, Hit]:
    """Return best Diamond hit per query transcript (max bitscore)."""
    best: dict[str, Hit] = {}
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
                evalue   = float(c[10])
                bitscore = float(c[11])
                qlen     = int(c[12])
            except (ValueError, IndexError):
                continue
            qcov = (qend - qstart + 1) / qlen * 100.0 if qlen > 0 else 0.0
            hit  = Hit(pident, qcov, evalue, bitscore)
            if qseqid not in best or bitscore > best[qseqid].bitscore:
                best[qseqid] = hit
    return best


def load_labels(labels_tsv: Path) -> dict[str, str]:
    with open(labels_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return {row["transcript_id"]: row["class_code"] for row in reader}


# ── Metrics ───────────────────────────────────────────────────────────────────

def pr_metrics(
    all_ids: set[str],
    labels:  dict[str, str],
    kept:    set[str],
) -> tuple[float, float, float]:
    """Return (precision, recall, F1) for the given kept set."""
    n_tp_total = sum(1 for t in all_ids if labels.get(t, "u") in TP_CLASS_CODES)
    n_tp_kept  = sum(1 for t in kept    if labels.get(t, "u") in TP_CLASS_CODES)
    n_fp_kept  = sum(1 for t in kept    if labels.get(t, "u") not in TP_CLASS_CODES)

    precision = n_tp_kept / (n_tp_kept + n_fp_kept) if (n_tp_kept + n_fp_kept) > 0 else 0.0
    recall    = n_tp_kept / n_tp_total if n_tp_total > 0 else 0.0
    f1        = (2 * precision * recall / (precision + recall)
                 if (precision + recall) > 0 else 0.0)
    return precision, recall, f1


# ── Threshold sweep ───────────────────────────────────────────────────────────

def sweep(
    all_ids: set[str],
    hits:    dict[str, Hit],
    labels:  dict[str, str],
) -> list[dict]:
    rows: list[dict] = []

    def record(rule: str, threshold, kept: set[str]) -> None:
        p, r, f1 = pr_metrics(all_ids, labels, kept)
        rows.append({
            "filter_rule": rule,
            "threshold":   threshold,
            "precision":   f"{p:.4f}",
            "recall":      f"{r:.4f}",
            "f1":          f"{f1:.4f}",
            "n_kept":      len(kept),
            "n_total":     len(all_ids),
        })

    # Baseline: no filtering
    record("no_filter", 0, all_ids)

    # Keep any transcript that has at least one hit
    has_hit = {t for t in all_ids if t in hits}
    record("has_any_hit", 0, has_hit)

    # Sweep pident (0–100, step 5)
    for thr in range(0, 101, 5):
        kept = {t for t, h in hits.items() if t in all_ids and h.pident >= thr}
        record("pident_ge", thr, kept)

    # Sweep query coverage (0–100, step 5)
    for thr in range(0, 101, 5):
        kept = {t for t, h in hits.items() if t in all_ids and h.query_cov >= thr}
        record("qcov_ge", thr, kept)

    # Sweep bitscore
    for thr in [0, 25, 50, 75, 100, 150, 200, 250, 300, 400, 500]:
        kept = {t for t, h in hits.items() if t in all_ids and h.bitscore >= thr}
        record("bitscore_ge", thr, kept)

    # Sweep e-value (report as -log10 for readability)
    for exp in [10, 5, 3, 2, 1, 0]:
        thr = 10 ** (-exp)
        kept = {t for t, h in hits.items() if t in all_ids and h.evalue <= thr}
        record("evalue_le_1e", -exp, kept)

    # Combined: pident ≥ X  AND  query_cov ≥ 50
    for pident_thr in range(0, 101, 10):
        kept = {
            t for t, h in hits.items()
            if t in all_ids and h.pident >= pident_thr and h.query_cov >= 50.0
        }
        record("pident_ge_qcov50", pident_thr, kept)

    # Combined: pident ≥ X  AND  query_cov ≥ 70
    for pident_thr in range(0, 101, 10):
        kept = {
            t for t, h in hits.items()
            if t in all_ids and h.pident >= pident_thr and h.query_cov >= 70.0
        }
        record("pident_ge_qcov70", pident_thr, kept)

    return rows


# ── Main ──────────────────────────────────────────────────────────────────────

def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--diamond",  required=True,
                    help="Diamond blastp TSV (14 cols, including qlen/slen)")
    ap.add_argument("--labels",   required=True,
                    help="Labels TSV from extract_labels.py")
    ap.add_argument("--species",  required=True)
    ap.add_argument("--level",    required=True,
                    help="Exclusion level label (genus / order / subphylum)")
    ap.add_argument("--out-dir",  required=True, type=Path)
    args = ap.parse_args(argv)

    args.out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[{args.species}/{args.level}] Loading Diamond hits ...", flush=True)
    hits = load_best_hits(Path(args.diamond))
    print(f"  -> {len(hits):,} queries with ≥1 hit", flush=True)

    print(f"[{args.species}/{args.level}] Loading labels ...", flush=True)
    labels  = load_labels(Path(args.labels))
    all_ids = set(labels)
    print(f"  -> {len(all_ids):,} transcripts with labels", flush=True)

    n_tp = sum(1 for t in all_ids if labels[t] in TP_CLASS_CODES)
    n_fp = len(all_ids) - n_tp
    print(f"  -> TP (codes {set(TP_CLASS_CODES)}): {n_tp:,}  |  FP: {n_fp:,}", flush=True)

    rows = sweep(all_ids, hits, labels)

    out_tsv = args.out_dir / f"{args.species}_{args.level}_pr_curves.tsv"
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
