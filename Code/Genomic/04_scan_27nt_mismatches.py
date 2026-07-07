#!/usr/bin/env python3
"""
Genome-wide mismatch scan for a 27-nt peptide-coding sequence (hg38)

Description
-----------
This script scans a reference genome FASTA (e.g., hg38 / GRCh38 primary assembly) for occurrences of a user-specified DNA query sequence (e.g., 27 nt) and reports the number of genomic matches with exactly 0..N mismatches (Hamming distance), considering both the query and its reverse complement.

It also writes BED6 files listing the genomic coordinates of hits for each mismatch class.

Output
------
- Console summary:
    * hits with exactly 0, 1, ..., N mismatches
    * hits with ≤0, ≤1, ..., ≤N mismatches
- BED6 files:
    <OUT_PREFIX>_mm0.bed
    <OUT_PREFIX>_mm1.bed
    ...
Coordinates are 0-based, end-exclusive (BED format).

Requirements
------------
- Python >= 3.8
- Reference genome FASTA (e.g., hg38.fa)

Example
-------
python3 scan_27mer_mismatches_hg38.py \
    --fasta hg38.fa \
    --query TTTGATCCTCAGCGGTCTCATTAGTAA \
    --max-mm 2 \
    --out-prefix LLM_27mer
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from typing import Dict, List, Tuple

Hit = Tuple[str, int, int, str]  # (chromosome, start0, end0, strand)


def revcomp(seq: str) -> str:
    """Return reverse complement of DNA sequence."""
    comp = str.maketrans("ACGTN", "TGCAN")
    return seq.translate(comp)[::-1]


def hamming_with_early_stop(a: str, b: str, max_mm: int) -> int:
    """Return Hamming distance between equal-length strings a and b with early stop."""
    mm = 0
    for x, y in zip(a, b):
        if x != y:
            mm += 1
            if mm > max_mm:
                return mm
    return mm


def flush_chromosome(
    chrom: str,
    seq_parts: List[str],
    query: str,
    query_rc: str,
    max_mm: int,
    hits_by_mm: Dict[int, List[Hit]],
) -> None:
    """Scan one chromosome sequence for hits up to max_mm mismatches."""
    seq = "".join(seq_parts).upper()
    k = len(query)
    n = len(seq)

    if n < k:
        return

    for i in range(0, n - k + 1):
        window = seq[i : i + k]

        if "N" in window:
            continue

        mm_plus = hamming_with_early_stop(window, query, max_mm)
        if mm_plus <= max_mm:
            hits_by_mm[mm_plus].append((chrom, i, i + k, "+"))

        mm_minus = hamming_with_early_stop(window, query_rc, max_mm)
        if mm_minus <= max_mm:
            hits_by_mm[mm_minus].append((chrom, i, i + k, "-"))


def write_bed6(path: str, hits: List[Hit], name: str) -> None:
    """Write hits in BED6 format."""
    with open(path, "w") as out:
        for chrom, start0, end0, strand in hits:
            out.write(f"{chrom}\t{start0}\t{end0}\t{name}\t0\t{strand}\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Genome-wide mismatch scan for a DNA query sequence."
    )
    parser.add_argument(
        "--fasta",
        required=True,
        help="Path to reference genome FASTA (e.g., hg38.fa).",
    )
    parser.add_argument(
        "--query",
        required=True,
        help="DNA query sequence (A/C/G/T only).",
    )
    parser.add_argument(
        "--max-mm",
        type=int,
        default=2,
        help="Maximum mismatches to report (default: 2).",
    )
    parser.add_argument(
        "--out-prefix",
        default="query",
        help="Prefix for output BED files (default: query).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    query = args.query.strip().upper()

    if any(base not in "ACGT" for base in query):
        raise SystemExit("ERROR: Query sequence must contain only A/C/G/T.")

    if args.max_mm < 0:
        raise SystemExit("ERROR: --max-mm must be >= 0.")

    query_rc = revcomp(query)
    max_mm = args.max_mm

    hits_by_mm: Dict[int, List[Hit]] = defaultdict(list)

    chrom = None
    seq_parts: List[str] = []

    with open(args.fasta, "r") as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if chrom is not None:
                    flush_chromosome(
                        chrom, seq_parts, query, query_rc, max_mm, hits_by_mm
                    )
                chrom = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)

        if chrom is not None:
            flush_chromosome(
                chrom, seq_parts, query, query_rc, max_mm, hits_by_mm
            )

    print(f"Query length: {len(query)}")
    print(f"Query (+): {query}")
    print(f"Query (-) revcomp: {query_rc}")
    print(f"Max mismatches scanned: {max_mm}")

    for mm in range(0, max_mm + 1):
        print(f"Hits with exactly {mm} mismatch(es): {len(hits_by_mm[mm])}")

    for mm in range(0, max_mm + 1):
        leq = sum(len(hits_by_mm[i]) for i in range(0, mm + 1))
        print(f"Hits with ≤{mm} mismatch(es): {leq}")

    for mm in range(0, max_mm + 1):
        out_bed = f"{args.out_prefix}_mm{mm}.bed"
        name = f"{args.out_prefix}_mm{mm}"
        write_bed6(out_bed, hits_by_mm[mm], name)
        print(f"Wrote {out_bed}")


if __name__ == "__main__":
    main()