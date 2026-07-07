#!/usr/bin/env bash
set -euo pipefail

# 01_EVA_mapping_to_hg38.sh
#
# Map EVA transcript sequences to the human reference genome (hg38) using
# splice-aware minimap2 alignment, then summarize transcript-level loci and
# exon-level genomic intervals from the alignment CIGAR strings.
#
# Intended use: Hulen et al. 2026 code availability/reproducibility repository.
#
# Requirements:
#   - bash
#   - python3 with PyYAML
#   - minimap2
#   - samtools
#   - bedtools
#
# Usage:
#   bash code/01_EVA_mapping_to_hg38.sh config/config.yaml
#
# Main outputs:
#   <outdir>/EVA_hg38.bam
#   <outdir>/EVA_hg38.bam.bai
#   <outdir>/hg38.genome
#   <outdir>/EVA_transcript_loci.tsv
#   <outdir>/EVA_exons_hg38.bed
#   <outdir>/EVA_exons_hg38.sorted.bed

CFG="${1:-config/config.yaml}"

if [[ ! -f "${CFG}" ]]; then
  echo "ERROR: Config file not found: ${CFG}" >&2
  exit 1
fi

read_yaml() {
  local key="$1"
  python3 - "$CFG" "$key" <<'PY'
import sys
import yaml

cfg_path, dotted_key = sys.argv[1], sys.argv[2]
with open(cfg_path, "r", encoding="utf-8") as handle:
    cfg = yaml.safe_load(handle)

value = cfg
for part in dotted_key.split("."):
    value = value[part]
print(value)
PY
}

EVA_FASTA="$(read_yaml inputs.eva_fasta)"
REF="$(read_yaml references.hg38_fasta)"
OUTDIR="$(read_yaml outputs.outdir)"
THREADS="$(read_yaml parameters.threads)"

for tool in minimap2 samtools bedtools python3; do
  if ! command -v "${tool}" >/dev/null 2>&1; then
    echo "ERROR: Required tool not found in PATH: ${tool}" >&2
    exit 1
  fi
done

if [[ ! -f "${EVA_FASTA}" ]]; then
  echo "ERROR: EVA FASTA not found: ${EVA_FASTA}" >&2
  exit 1
fi

if [[ ! -f "${REF}" ]]; then
  echo "ERROR: hg38 reference FASTA not found: ${REF}" >&2
  exit 1
fi

mkdir -p "${OUTDIR}"

# Ensure reference is indexed for downstream genomic scanning and contig ordering.
if [[ ! -f "${REF}.fai" ]]; then
  samtools faidx "${REF}"
fi

cut -f1,2 "${REF}.fai" > "${OUTDIR}/hg38.genome"

# Splice-aware transcript-to-genome alignment.
minimap2 -t "${THREADS}" -ax splice -uf -k14 "${REF}" "${EVA_FASTA}" \
  | samtools sort -@ "${THREADS}" -o "${OUTDIR}/EVA_hg38.bam"

samtools index "${OUTDIR}/EVA_hg38.bam"

# Summarize transcript loci and exon blocks. This intentionally reads SAM records
# via subprocess rather than piping into a heredoc, so stdin remains unambiguous.
python3 - "${OUTDIR}/EVA_hg38.bam" "${OUTDIR}/EVA_transcript_loci.tsv" "${OUTDIR}/EVA_exons_hg38.bed" <<'PY'
import re
import subprocess
import sys
from pathlib import Path

bam_path = Path(sys.argv[1])
loci_out = Path(sys.argv[2])
exon_out = Path(sys.argv[3])

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def cigar_ref_span(cigar: str) -> int:
    total = 0
    for n, op in CIGAR_RE.findall(cigar):
        if op in {"M", "D", "N", "=", "X"}:
            total += int(n)
    return total


def cigar_exons(pos_1based: int, cigar: str):
    """Return 0-based half-open exon intervals split at skipped regions (N)."""
    ref_pos_0 = pos_1based - 1
    exon_start = ref_pos_0
    exon_end = ref_pos_0
    exons = []

    for n_str, op in CIGAR_RE.findall(cigar):
        n = int(n_str)
        if op in {"M", "=", "X", "D"}:
            exon_end += n
            ref_pos_0 += n
        elif op == "N":
            if exon_end > exon_start:
                exons.append((exon_start, exon_end))
            ref_pos_0 += n
            exon_start = ref_pos_0
            exon_end = ref_pos_0
        elif op in {"I", "S", "H", "P"}:
            # Does not consume reference. Soft clipping consumes query only.
            continue

    if exon_end > exon_start:
        exons.append((exon_start, exon_end))
    return exons


sam = subprocess.run(
    ["samtools", "view", str(bam_path)],
    check=True,
    text=True,
    capture_output=True,
).stdout.splitlines()

with loci_out.open("w", encoding="utf-8") as loci, exon_out.open("w", encoding="utf-8") as exons_handle:
    loci.write("Transcript\tChrom\tStart_1based\tEnd_1based\tStrand\tMAPQ\tCIGAR\n")
    for line in sam:
        fields = line.rstrip("\n").split("\t")
        qname = fields[0]
        flag = int(fields[1])
        chrom = fields[2]
        pos = int(fields[3])
        mapq = fields[4]
        cigar = fields[5]

        if chrom == "*" or cigar == "*":
            continue

        strand = "-" if (flag & 16) else "+"
        end = pos + cigar_ref_span(cigar) - 1
        loci.write(f"{qname}\t{chrom}\t{pos}\t{end}\t{strand}\t{mapq}\t{cigar}\n")

        for idx, (start0, end0) in enumerate(cigar_exons(pos, cigar), start=1):
            exons_handle.write(f"{chrom}\t{start0}\t{end0}\t{qname}_exon{idx}\t{mapq}\t{strand}\n")
PY

bedtools sort -g "${OUTDIR}/hg38.genome" -i "${OUTDIR}/EVA_exons_hg38.bed" \
  > "${OUTDIR}/EVA_exons_hg38.sorted.bed"

echo "Completed EVA mapping and exon extraction. Outputs written to: ${OUTDIR}"
