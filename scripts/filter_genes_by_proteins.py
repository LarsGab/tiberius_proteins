#!/usr/bin/env python3
"""Adaptive protein-based filtering of a predicted gene set.

Pipeline:
    1. Extract peptide sequences from <genome> + <gtf>          (gffread)
    2. Build Diamond DB from <protein FASTA>                    (diamond makedb)
    3. Run Diamond blastp                                       (diamond blastp)
    4. Estimate DB-to-genome distance from the best-hit %ID
       distribution and the fraction of queries with any hit.
    5. Pick a FilterRule (default: AdaptiveRule from filter_rules.py).
    6. Apply rule -> keep only transcripts with a passing hit
       (or all transcripts if the rule is disabled).
    7. Write filtered GTF + filter_report.json.

Existing intermediate files are reused unless --force is given.

Usage:
    python filter_genes_by_proteins.py \\
        --genome    genome.fa \\
        --gtf       predictions.gtf \\
        --proteins  protein_db.fa \\
        --out-dir   results/filtered/ \\
        [--threads 8]

Override the adaptive rule with a fixed one:
    --rule fixed --pident-min 50 --qcov-min 50 --evalue-max 1e-5

Skip the expensive Diamond step if you already have its output:
    --diamond-tsv path/to/diamond.tsv

To change the adaptive defaults (close_threshold, per-bin rules, …)
edit scripts/filter_rules.py.
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import statistics
import subprocess
import sys
from dataclasses import asdict
from pathlib import Path
from typing import Iterable, Optional

# Make filter_rules importable when this script is invoked from anywhere
sys.path.insert(0, str(Path(__file__).resolve().parent))
from filter_rules import (
    AdaptiveRule,
    DBStats,
    FilterRule,
    Hit,
)

# Diamond blastp output columns we request
_DIAMOND_FMT = (
    "qseqid sseqid pident length mismatch gapopen "
    "qstart qend sstart send evalue bitscore qlen slen"
)
_TRANSCRIPT_ID_RE = re.compile(r'transcript_id "([^"]+)"')
_GENE_ID_RE       = re.compile(r'gene_id "([^"]+)"')


# ── Subprocess wrappers ───────────────────────────────────────────────────────
def _run(cmd: list[str], **kwargs) -> None:
    print(f"  $ {' '.join(cmd)}", flush=True)
    subprocess.run(cmd, check=True, **kwargs)


def extract_peptides(
    genome: Path,
    gtf: Path,
    out_fa: Path,
    *,
    force: bool = False,
) -> Path:
    """Translate CDS in `gtf` to protein sequences via gffread."""
    if out_fa.exists() and not force:
        print(f"[1/7] peptides exist, reusing: {out_fa}", flush=True)
        return out_fa
    if shutil.which("gffread") is None:
        raise RuntimeError("gffread not found in PATH")

    print(f"[1/7] extracting peptides from {gtf}", flush=True)
    out_fa.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_fa.with_suffix(".tmp.fa")
    _run(["gffread", "-y", str(tmp), "-g", str(genome), str(gtf)])

    # strip stop-codon asterisks and gffread masked-residue dots; Diamond rejects both
    with open(tmp) as fin, open(out_fa, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                fout.write(line)
            else:
                fout.write(line.replace("*", "").replace(".", ""))
    tmp.unlink()
    return out_fa


def run_diamond(
    query_fa: Path,
    proteins_fa: Path,
    out_dir: Path,
    *,
    threads: int = 8,
    evalue: float = 1e-5,
    sensitivity: str = "--more-sensitive",
    max_target_seqs: int = 5,
    force: bool = False,
) -> Path:
    """Run `diamond makedb` then `blastp`. Returns the path of the blastp TSV."""
    if shutil.which("diamond") is None:
        raise RuntimeError("diamond not found in PATH")

    db_prefix = out_dir / "diamond_db"
    out_tsv   = out_dir / "diamond.tsv"
    out_dir.mkdir(parents=True, exist_ok=True)

    if out_tsv.exists() and not force:
        print(f"[2-3/7] diamond.tsv exists, reusing: {out_tsv}", flush=True)
        return out_tsv

    if not db_prefix.with_suffix(".dmnd").exists() or force:
        print(f"[2/7] building Diamond DB", flush=True)
        _run([
            "diamond", "makedb",
            "--in",     str(proteins_fa),
            "--db",     str(db_prefix),
            "--threads", str(threads),
        ])

    print(f"[3/7] running diamond blastp", flush=True)
    _run([
        "diamond", "blastp",
        "--query",   str(query_fa),
        "--db",      str(db_prefix.with_suffix(".dmnd")),
        "--out",     str(out_tsv),
        "--outfmt",  "6", *_DIAMOND_FMT.split(),
        "--max-target-seqs", str(max_target_seqs),
        "--evalue",  str(evalue),
        sensitivity,
        "--threads", str(threads),
    ])
    return out_tsv


# ── Diamond TSV parsing ───────────────────────────────────────────────────────
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
                length   = int(c[3])
                qstart   = int(c[6])
                qend     = int(c[7])
                evalue   = float(c[10])
                bitscore = float(c[11])
                qlen     = int(c[12])
                slen     = int(c[13])
            except (ValueError, IndexError):
                continue
            qcov = (qend - qstart + 1) / qlen * 100.0 if qlen > 0 else 0.0
            tcov = length / slen * 100.0 if slen > 0 else 0.0
            hit = Hit(pident, qcov, tcov, evalue, bitscore)
            if qseqid not in best or bitscore > best[qseqid].bitscore:
                best[qseqid] = hit
    return best


# ── GTF helpers ───────────────────────────────────────────────────────────────
def gtf_transcript_ids(gtf: Path) -> set[str]:
    ids: set[str] = set()
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            m = _TRANSCRIPT_ID_RE.search(line)
            if m:
                ids.add(m.group(1))
    return ids


def write_filtered_gtf(
    src_gtf: Path,
    dst_gtf: Path,
    kept_transcript_ids: set[str],
) -> tuple[int, int]:
    """Write `src_gtf` lines whose transcript_id is in `kept_transcript_ids`.

    Gene-level lines (no transcript_id) are kept iff their gene_id has at
    least one kept transcript.
    Returns (n_lines_written, n_lines_total).
    """
    # Pass 1: which gene_ids are kept?
    kept_gene_ids: set[str] = set()
    with open(src_gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            tid_m = _TRANSCRIPT_ID_RE.search(line)
            gid_m = _GENE_ID_RE.search(line)
            if tid_m and gid_m and tid_m.group(1) in kept_transcript_ids:
                kept_gene_ids.add(gid_m.group(1))

    # Pass 2: write filtered output
    n_kept = n_total = 0
    with open(src_gtf) as fin, open(dst_gtf, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            n_total += 1
            tid_m = _TRANSCRIPT_ID_RE.search(line)
            gid_m = _GENE_ID_RE.search(line)
            keep = False
            if tid_m:
                keep = tid_m.group(1) in kept_transcript_ids
            elif gid_m:
                keep = gid_m.group(1) in kept_gene_ids
            if keep:
                fout.write(line)
                n_kept += 1
    return n_kept, n_total


# ── Distance estimation ───────────────────────────────────────────────────────
def estimate_db_stats(query_ids: Iterable[str], hits: dict[str, Hit]) -> DBStats:
    query_ids = set(query_ids)
    hit_pidents = [h.pident for tid, h in hits.items() if tid in query_ids]
    return DBStats(
        n_queries     = len(query_ids),
        n_with_hit    = len(hit_pidents),
        median_pident = statistics.median(hit_pidents) if hit_pidents else 0.0,
    )


# ── Rule construction from CLI ────────────────────────────────────────────────
def build_fixed_rule(args: argparse.Namespace) -> FilterRule:
    return FilterRule(
        name         = "cli_fixed",
        pident_min   = args.pident_min,
        qcov_min     = args.qcov_min,
        tcov_min     = args.tcov_min,
        evalue_max   = args.evalue_max,
        bitscore_min = args.bitscore_min,
        disabled     = False,
    )


# ── Main ──────────────────────────────────────────────────────────────────────
def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    inp = ap.add_argument_group("inputs")
    inp.add_argument("--genome",   type=Path, required=True)
    inp.add_argument("--gtf",      type=Path, required=True)
    inp.add_argument("--proteins", type=Path, required=True)

    out = ap.add_argument_group("outputs")
    out.add_argument("--out-dir",  type=Path, required=True)

    cache = ap.add_argument_group("caching / reuse")
    cache.add_argument("--peptides-fa", type=Path,
                       help="Skip gffread, use this peptide FASTA")
    cache.add_argument("--diamond-tsv", type=Path,
                       help="Skip Diamond, use this existing TSV (must be 14-col format)")
    cache.add_argument("--force", action="store_true",
                       help="Recompute every step even if outputs already exist")

    diamond = ap.add_argument_group("diamond")
    diamond.add_argument("--threads", type=int, default=8)
    diamond.add_argument("--diamond-evalue", type=float, default=1e-5)
    diamond.add_argument("--sensitivity",
                         choices=["--fast", "--mid-sensitive", "--sensitive",
                                  "--more-sensitive", "--very-sensitive", "--ultra-sensitive"],
                         default="--more-sensitive")
    diamond.add_argument("--max-target-seqs", type=int, default=5)

    rule = ap.add_argument_group("filter rule")
    rule.add_argument("--rule", choices=["adaptive", "fixed"], default="adaptive",
                      help="adaptive (default): pick rule by median pident; "
                           "fixed: use the thresholds below")
    rule.add_argument("--pident-min",   type=float, default=0.0)
    rule.add_argument("--qcov-min",     type=float, default=0.0)
    rule.add_argument("--tcov-min",     type=float, default=0.0)
    rule.add_argument("--evalue-max",   type=float, default=float("inf"))
    rule.add_argument("--bitscore-min", type=float, default=0.0)

    return ap.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: peptides
    peptides_fa = args.peptides_fa
    if peptides_fa is None:
        peptides_fa = extract_peptides(
            args.genome, args.gtf, args.out_dir / "peptides.fa", force=args.force
        )

    # Steps 2–3: Diamond
    diamond_tsv = args.diamond_tsv
    if diamond_tsv is None:
        diamond_tsv = run_diamond(
            peptides_fa, args.proteins, args.out_dir,
            threads=args.threads,
            evalue=args.diamond_evalue,
            sensitivity=args.sensitivity,
            max_target_seqs=args.max_target_seqs,
            force=args.force,
        )

    # Step 4: best hits + DB stats
    print(f"[4/7] loading best hits", flush=True)
    best_hits = load_best_hits(diamond_tsv)
    all_tids  = gtf_transcript_ids(args.gtf)
    stats     = estimate_db_stats(all_tids, best_hits)
    print(
        f"  queries={stats.n_queries:,}  with_hit={stats.n_with_hit:,} "
        f"({stats.hit_fraction:.1%})  median best-hit pident={stats.median_pident:.2f}%",
        flush=True,
    )

    # Step 5: pick rule
    print(f"[5/7] selecting filter rule", flush=True)
    if args.rule == "fixed":
        rule = build_fixed_rule(args)
    else:
        rule = AdaptiveRule().select(stats)
    print(f"  -> {rule.describe()}", flush=True)

    # Step 6: apply rule
    print(f"[6/7] applying rule to {stats.n_queries:,} transcripts", flush=True)
    kept_tids: set[str] = set()
    for tid in all_tids:
        if rule.passes(best_hits.get(tid)):
            kept_tids.add(tid)
    print(
        f"  kept {len(kept_tids):,} / {stats.n_queries:,} "
        f"({len(kept_tids) / max(1, stats.n_queries):.1%})",
        flush=True,
    )

    # Step 7: write outputs
    out_gtf = args.out_dir / "filtered.gtf"
    n_kept_lines, n_total_lines = write_filtered_gtf(args.gtf, out_gtf, kept_tids)
    print(f"[7/7] wrote {out_gtf}  ({n_kept_lines:,} / {n_total_lines:,} lines)",
          flush=True)

    report = {
        "inputs": {
            "genome":   str(args.genome),
            "gtf":      str(args.gtf),
            "proteins": str(args.proteins),
        },
        "diamond_tsv":   str(diamond_tsv),
        "db_stats":      asdict(stats),
        "rule": {
            "type":        args.rule,
            "name":        rule.name,
            "describe":    rule.describe(),
            "pident_min":  rule.pident_min,
            "qcov_min":    rule.qcov_min,
            "tcov_min":    rule.tcov_min,
            "evalue_max":  rule.evalue_max,
            "bitscore_min":rule.bitscore_min,
            "disabled":    rule.disabled,
        },
        "result": {
            "n_transcripts_total": stats.n_queries,
            "n_transcripts_kept":  len(kept_tids),
            "n_lines_total":       n_total_lines,
            "n_lines_kept":        n_kept_lines,
        },
    }
    report_path = args.out_dir / "filter_report.json"
    report_path.write_text(json.dumps(report, indent=2, default=str))
    print(f"      report : {report_path}", flush=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
