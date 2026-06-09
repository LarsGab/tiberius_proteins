#!/usr/bin/env python3
"""Locus-style Diamond filter for Tiberius predictions.

Replicates the miniprot+miniprothint filter (keep a Tiberius gene if its
locus has any protein hint) with Diamond hits instead: a gene is kept iff
at least one of its transcripts has a Diamond hit above a minimum e-value.
No pident / qcov tuning -- it is just "has any hit".

Two input modes:
  (A) provide a protein DB FASTA (and genome) -> the script runs gffread
      + diamond makedb + diamond blastp itself
  (B) provide a Diamond blastp TSV (outfmt 6, at least qseqid + evalue)
      -> Diamond step is skipped

The query IDs in the Diamond TSV must match the transcript_id values in
the GTF (gffread emits transcript_id as the FASTA header, which is what
the existing pipeline uses).

Usage (mode A):
    filter_tiberius_by_diamond.py \\
        --gtf       tiberius.gtf \\
        --genome    genome.fa \\
        --proteins  protein_db.fa \\
        --out-dir   results/diamond_filter/ \\
        [--threads 8] [--evalue 1e-5]

Usage (mode B):
    filter_tiberius_by_diamond.py \\
        --gtf         tiberius.gtf \\
        --diamond-tsv diamond.tsv \\
        --out-dir     results/diamond_filter/ \\
        [--evalue 1e-5]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

# Reuse helpers from the existing transcript-level filter script
sys.path.insert(0, str(Path(__file__).resolve().parent))
from filter_genes_by_proteins import (  # noqa: E402
    _GENE_ID_RE,
    _TRANSCRIPT_ID_RE,
    extract_peptides,
    run_diamond,
    write_filtered_gtf,
)


from dataclasses import dataclass

_SCI_NOTATION_HINT = ("e", "E")


def _looks_like_evalue(tok: str) -> bool:
    """A token is e-value-like if it parses as float AND uses sci notation."""
    if not any(c in tok for c in _SCI_NOTATION_HINT):
        return False
    try:
        float(tok)
    except ValueError:
        return False
    return True


def detect_evalue_col(diamond_tsv: Path, scan_lines: int = 50) -> int:
    """Return the 0-based column index that looks like the e-value column.

    Heuristic: find the column whose tokens consistently use scientific
    notation (e.g. '1e-5', '3.33e-239'). Diamond writes bitscore as a plain
    number, so this uniquely identifies e-value across common outfmt 6
    variants (default 12-col -> idx 10; custom layouts -> different idx).
    """
    counts: dict[int, int] = {}
    n_lines = 0
    with open(diamond_tsv) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            for i, tok in enumerate(cols):
                if _looks_like_evalue(tok):
                    counts[i] = counts.get(i, 0) + 1
            n_lines += 1
            if n_lines >= scan_lines:
                break
    if not counts:
        raise ValueError(
            f"Could not auto-detect e-value column in {diamond_tsv}. "
            "Pass --evalue-col explicitly (1-based)."
        )
    # Column with the most sci-notation tokens wins; tiebreak on lowest index.
    best_idx = max(counts, key=lambda i: (counts[i], -i))
    return best_idx


@dataclass
class ColLayout:
    """0-based column indices in the Diamond TSV.

    Both common Diamond outfmt 6 variants are handled with the same offsets
    from the e-value column:
        14-col extended : pident=2 length=3  evalue=10 bitscore=11 qlen=12 slen=13
        8-col compact   : pident=2 length=3  evalue=4  bitscore=5  qlen=6  slen=7
    bitscore = evalue+1, qlen = evalue+2, slen = evalue+3 in both layouts.
    """
    pident:   int
    length:   int
    evalue:   int
    bitscore: int
    qlen:     int  # -1 if not present
    slen:     int  # -1 if not present

    @classmethod
    def from_evalue_col(cls, evalue_col: int, n_cols: int) -> "ColLayout":
        bitscore = evalue_col + 1 if evalue_col + 1 < n_cols else -1
        qlen     = evalue_col + 2 if evalue_col + 2 < n_cols else -1
        slen     = evalue_col + 3 if evalue_col + 3 < n_cols else -1
        return cls(pident=2, length=3, evalue=evalue_col,
                   bitscore=bitscore, qlen=qlen, slen=slen)


@dataclass
class HitFilter:
    """Per-hit thresholds. A hit passes iff ALL active thresholds are met."""
    pident_min:   float = 0.0
    qcov_min:     float = 0.0
    tcov_min:     float = 0.0
    evalue_max:   float = float("inf")
    bitscore_min: float = 0.0
    length_min:   int   = 0

    def describe(self) -> str:
        p = []
        if self.pident_min   > 0:            p.append(f"pident>={self.pident_min:g}")
        if self.qcov_min     > 0:            p.append(f"qcov>={self.qcov_min:g}")
        if self.tcov_min     > 0:            p.append(f"tcov>={self.tcov_min:g}")
        if self.evalue_max   < float("inf"): p.append(f"evalue<={self.evalue_max:g}")
        if self.bitscore_min > 0:            p.append(f"bitscore>={self.bitscore_min:g}")
        if self.length_min   > 0:            p.append(f"length>={self.length_min}")
        return ", ".join(p) if p else "any-hit"


def transcripts_with_hit(
    diamond_tsv: Path,
    flt: HitFilter,
    layout: ColLayout,
) -> tuple[set[str], int, int]:
    """Return (supported_qids, n_rows_seen, n_rows_passing)."""
    keep: set[str] = set()
    n_rows = n_pass = 0
    with open(diamond_tsv) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) <= layout.evalue:
                continue
            n_rows += 1
            try:
                pident = float(cols[layout.pident])
                length = int(cols[layout.length])
                evalue = float(cols[layout.evalue])
                bitscore = (float(cols[layout.bitscore])
                            if layout.bitscore >= 0 else 0.0)
            except (ValueError, IndexError):
                continue

            if pident   < flt.pident_min:   continue
            if evalue   > flt.evalue_max:   continue
            if bitscore < flt.bitscore_min: continue
            if length   < flt.length_min:   continue

            if flt.qcov_min > 0:
                if layout.qlen < 0:
                    sys.exit("ERROR: --qcov-min set but qlen column not "
                             "present in TSV. Run Diamond with qlen in "
                             "--outfmt, or drop --qcov-min.")
                try:
                    qlen = int(cols[layout.qlen])
                except (ValueError, IndexError):
                    continue
                qcov = length / qlen * 100.0 if qlen > 0 else 0.0
                if qcov < flt.qcov_min:
                    continue

            if flt.tcov_min > 0:
                if layout.slen < 0:
                    sys.exit("ERROR: --tcov-min set but slen column not "
                             "present in TSV.")
                try:
                    slen = int(cols[layout.slen])
                except (ValueError, IndexError):
                    continue
                tcov = length / slen * 100.0 if slen > 0 else 0.0
                if tcov < flt.tcov_min:
                    continue

            keep.add(cols[0])
            n_pass += 1
    return keep, n_rows, n_pass


def gene_to_transcripts(gtf: Path) -> dict[str, set[str]]:
    """Map gene_id -> set of transcript_ids found in the GTF."""
    mapping: dict[str, set[str]] = {}
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            tid = _TRANSCRIPT_ID_RE.search(line)
            gid = _GENE_ID_RE.search(line)
            if tid and gid:
                mapping.setdefault(gid.group(1), set()).add(tid.group(1))
    return mapping


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--gtf", required=True, type=Path,
                    help="Tiberius predictions (GTF)")
    ap.add_argument("--out-dir", required=True, type=Path)

    # Mode A: run Diamond ourselves
    ap.add_argument("--genome", type=Path,
                    help="Genome FASTA (required if --diamond-tsv not given)")
    ap.add_argument("--proteins", type=Path,
                    help="Protein DB FASTA (required if --diamond-tsv not given)")
    ap.add_argument("--peptides-fa", type=Path,
                    help="Skip gffread, use this peptide FASTA as Diamond query")

    # Mode B: pre-computed Diamond TSV
    ap.add_argument("--diamond-tsv", type=Path,
                    help="Pre-computed Diamond blastp TSV (outfmt 6). "
                         "If given, --genome/--proteins are ignored.")

    # Filtering thresholds
    flt = ap.add_argument_group("filter thresholds (a hit counts as 'support' "
                                "iff ALL active thresholds are met)")
    flt.add_argument("--evalue",       type=float, default=1e-5,
                     help="Max e-value (default 1e-5)")
    flt.add_argument("--pident-min",   type=float, default=0.0,
                     help="Min %% identity, 0-100 (default 0 = off)")
    flt.add_argument("--qcov-min",     type=float, default=0.0,
                     help="Min query coverage in %% (length/qlen*100). "
                          "Requires qlen column in TSV. Default 0 = off.")
    flt.add_argument("--tcov-min",     type=float, default=0.0,
                     help="Min target coverage in %% (length/slen*100). "
                          "Requires slen column in TSV. Default 0 = off.")
    flt.add_argument("--bitscore-min", type=float, default=0.0,
                     help="Min bitscore (default 0 = off)")
    flt.add_argument("--length-min",   type=int,   default=0,
                     help="Min alignment length (default 0 = off)")

    # Column overrides (auto-detected for evalue; bitscore/qlen/slen are
    # inferred as evalue+1/+2/+3)
    cols = ap.add_argument_group("column overrides (1-based)")
    cols.add_argument("--evalue-col",   type=int, default=None,
                      help="If omitted, auto-detected by scanning for sci "
                           "notation.")
    cols.add_argument("--pident-col",   type=int, default=3,
                      help="Default 3 (Diamond convention).")
    cols.add_argument("--length-col",   type=int, default=4,
                      help="Default 4 (Diamond convention).")
    cols.add_argument("--bitscore-col", type=int, default=None,
                      help="Default: evalue-col + 1.")
    cols.add_argument("--qlen-col",     type=int, default=None,
                      help="Default: evalue-col + 2 if present.")
    cols.add_argument("--slen-col",     type=int, default=None,
                      help="Default: evalue-col + 3 if present.")

    # Diamond knobs (only used in mode A)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--sensitivity",
                    choices=["--fast", "--mid-sensitive", "--sensitive",
                             "--more-sensitive", "--very-sensitive",
                             "--ultra-sensitive"],
                    default="--more-sensitive")
    ap.add_argument("--max-target-seqs", type=int, default=1,
                    help="Only need one hit per query (default 1)")
    ap.add_argument("--force", action="store_true",
                    help="Recompute peptides / Diamond even if outputs exist")

    return ap.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    if args.diamond_tsv is None:
        if args.genome is None or args.proteins is None:
            sys.exit("ERROR: provide either --diamond-tsv OR --genome + --proteins")

        peptides_fa = args.peptides_fa or extract_peptides(
            args.genome, args.gtf, args.out_dir / "peptides.fa", force=args.force,
        )
        diamond_tsv = run_diamond(
            peptides_fa, args.proteins, args.out_dir,
            threads=args.threads,
            evalue=args.evalue,
            sensitivity=args.sensitivity,
            max_target_seqs=args.max_target_seqs,
            force=args.force,
        )
    else:
        diamond_tsv = args.diamond_tsv
        print(f"Using pre-computed Diamond TSV: {diamond_tsv}", flush=True)

    # Resolve columns (CLI is 1-based; layout is 0-based)
    if args.evalue_col is None:
        evalue_col = detect_evalue_col(diamond_tsv)
        print(f"Auto-detected e-value column: {evalue_col + 1} (1-based)",
              flush=True)
    else:
        evalue_col = args.evalue_col - 1

    # Peek at first data row to know how many columns we have
    n_cols = 0
    with open(diamond_tsv) as fh:
        for line in fh:
            if line.strip() and not line.startswith("#"):
                n_cols = len(line.rstrip("\n").split("\t"))
                break
    layout = ColLayout.from_evalue_col(evalue_col, n_cols)
    layout.pident = args.pident_col - 1
    layout.length = args.length_col - 1
    if args.bitscore_col is not None: layout.bitscore = args.bitscore_col - 1
    if args.qlen_col     is not None: layout.qlen     = args.qlen_col - 1
    if args.slen_col     is not None: layout.slen     = args.slen_col - 1

    flt = HitFilter(
        pident_min   = args.pident_min,
        qcov_min     = args.qcov_min,
        tcov_min     = args.tcov_min,
        evalue_max   = args.evalue,
        bitscore_min = args.bitscore_min,
        length_min   = args.length_min,
    )
    print(f"Filter: {flt.describe()}", flush=True)
    supported_tids, n_rows, n_pass = transcripts_with_hit(
        diamond_tsv, flt, layout,
    )
    print(f"  rows scanned : {n_rows:,}", flush=True)
    print(f"  rows passing : {n_pass:,}", flush=True)

    gene_tids = gene_to_transcripts(args.gtf)
    n_genes_total = len(gene_tids)
    n_tids_total  = sum(len(v) for v in gene_tids.values())

    # Locus-style: a gene is kept if ANY of its transcripts has a hit.
    # All transcripts of a kept gene are then retained.
    kept_genes = {gid for gid, tids in gene_tids.items()
                  if tids & supported_tids}
    kept_tids: set[str] = set()
    for gid in kept_genes:
        kept_tids.update(gene_tids[gid])

    print(f"  transcripts with hit : {len(supported_tids):,} / {n_tids_total:,}",
          flush=True)
    print(f"  genes kept           : {len(kept_genes):,} / {n_genes_total:,} "
          f"({len(kept_genes)/max(1,n_genes_total):.1%})", flush=True)
    print(f"  transcripts kept     : {len(kept_tids):,} / {n_tids_total:,}",
          flush=True)

    out_gtf = args.out_dir / "filtered.gtf"
    n_kept_lines, n_total_lines = write_filtered_gtf(args.gtf, out_gtf, kept_tids)
    print(f"Wrote {out_gtf}  ({n_kept_lines:,} / {n_total_lines:,} lines)",
          flush=True)

    report = {
        "inputs": {
            "gtf":         str(args.gtf),
            "genome":      str(args.genome)   if args.genome   else None,
            "proteins":    str(args.proteins) if args.proteins else None,
            "diamond_tsv": str(diamond_tsv),
        },
        "filter": {
            "describe":     flt.describe(),
            "evalue_max":   flt.evalue_max,
            "pident_min":   flt.pident_min,
            "qcov_min":     flt.qcov_min,
            "tcov_min":     flt.tcov_min,
            "bitscore_min": flt.bitscore_min,
            "length_min":   flt.length_min,
        },
        "columns_1based": {
            "pident":   layout.pident + 1,
            "length":   layout.length + 1,
            "evalue":   layout.evalue + 1,
            "bitscore": layout.bitscore + 1 if layout.bitscore >= 0 else None,
            "qlen":     layout.qlen + 1     if layout.qlen     >= 0 else None,
            "slen":     layout.slen + 1     if layout.slen     >= 0 else None,
        },
        "diamond_rows_scanned": n_rows,
        "diamond_rows_passing": n_pass,
        "result": {
            "n_genes_total":         n_genes_total,
            "n_genes_kept":          len(kept_genes),
            "n_transcripts_total":   n_tids_total,
            "n_transcripts_with_hit":len(supported_tids),
            "n_transcripts_kept":    len(kept_tids),
            "n_lines_total":         n_total_lines,
            "n_lines_kept":          n_kept_lines,
        },
    }
    (args.out_dir / "filter_report.json").write_text(
        json.dumps(report, indent=2, default=str)
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
