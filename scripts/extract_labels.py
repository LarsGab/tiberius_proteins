#!/usr/bin/env python3
"""Extract per-transcript gffcompare class codes from an annotated GTF.

gffcompare adds a `class_code` attribute to every transcript line in the
annotated GTF it produces.  Common codes used here:

    =   exact match of all intron/exon boundaries (strict TP)
    c   contained within a reference transcript (intron-compatible)
    j   potential novel isoform (≥1 junction matches)
    u   unknown / intergenic (clear FP)
    i   fully intronic
    x   exonic overlap on opposite strand
    s   intron overlap on opposite strand
    e   single-exon transfrag partially covering a reference intron
    o   other exonic overlap

Output: two-column TSV  transcript_id <TAB> class_code

Usage:
    extract_labels.py ANNOTATED_GTF OUTPUT_TSV
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

_TRANSCRIPT_ID = re.compile(r'transcript_id\s+"([^"]+)"')
_CLASS_CODE    = re.compile(r'class_code\s+"([^"]+)"')


def parse_annotated_gtf(gtf_path: Path) -> dict[str, str]:
    labels: dict[str, str] = {}
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            attrs = fields[8]
            tid_m = _TRANSCRIPT_ID.search(attrs)
            cc_m  = _CLASS_CODE.search(attrs)
            if tid_m and cc_m:
                labels[tid_m.group(1)] = cc_m.group(1)
    return labels


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("annotated_gtf", help="gffcompare *.annotated.gtf")
    ap.add_argument("output_tsv",    help="Output TSV (transcript_id, class_code)")
    args = ap.parse_args(argv)

    labels = parse_annotated_gtf(Path(args.annotated_gtf))
    with open(args.output_tsv, "w") as fh:
        fh.write("transcript_id\tclass_code\n")
        for tid, cc in sorted(labels.items()):
            fh.write(f"{tid}\t{cc}\n")

    print(f"Wrote {len(labels):,} transcript labels to {args.output_tsv}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
