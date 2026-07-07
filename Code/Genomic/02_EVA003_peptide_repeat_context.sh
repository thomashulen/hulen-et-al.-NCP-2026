#!/usr/bin/env bash
set -euo pipefail

# 02_EVA003_peptide_repeat_context.sh
#
# Locate an epitope-coding peptide sequence within a target EVA transcript,
# project the transcript nucleotide interval to hg38 genomic coordinates using
# the spliced alignment from Script 01, and quantify RepeatMasker coverage over
# the target transcript locus.
#
# Intended use: Hulen et al. 2026 code availability/reproducibility repository.
#
# Requirements:
#   - bash
#   - python3 with PyYAML
#   - samtools
#   - bedtools
#
# Usage:
#   bash code/02_EVA003_peptide_repeat_context.sh config/config.yaml
#
# Requires Script 01 outputs in <outdir>:
#   EVA_hg38.bam
#   EVA_transcript_loci.tsv
#   hg38.genome
#
# Main outputs:
#   <outdir>/peptide_in_<target>_transcript.tsv
#   <outdir>/peptide_<peptide>_genomic_coords.tsv
#   <outdir>/<target>_locus_hg38.sorted.bed
#   <outdir>/<target>_repeat_composition_clean.tsv
#   <outdir>/<target>_locus_LTR_hits.tsv
#   <outdir>/<peptide>_27nt.fa

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
RMSK_BED="$(read_yaml references.rmsk_bed)"
TARGET_TX="$(read_yaml targets.transcript_id)"
PEP="$(read_yaml targets.peptide)"
OUTDIR="$(read_yaml outputs.outdir)"
RMSK_CLASS_COL="$(read_yaml parameters.rmsk_class_column)"

for tool in samtools bedtools python3; do
  if ! command -v "${tool}" >/dev/null 2>&1; then
    echo "ERROR: Required tool not found in PATH: ${tool}" >&2
    exit 1
  fi
done

for f in "${EVA_FASTA}" "${RMSK_BED}" "${OUTDIR}/EVA_hg38.bam" "${OUTDIR}/EVA_transcript_loci.tsv" "${OUTDIR}/hg38.genome"; do
  if [[ ! -f "${f}" ]]; then
    echo "ERROR: Required input not found: ${f}" >&2
    exit 1
  fi
done

# 1) Find peptide sequence within the target transcript in all three forward frames.
python3 - "${EVA_FASTA}" "${TARGET_TX}" "${PEP}" "${OUTDIR}/peptide_in_${TARGET_TX}_transcript.tsv" <<'PY'
from pathlib import Path
import sys

fasta = Path(sys.argv[1])
target = sys.argv[2]
peptide = sys.argv[3]
out_path = Path(sys.argv[4])

CODON_TABLE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K","GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}


def read_fasta_record(path: Path, record_id: str) -> str:
    seq = []
    found = False
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                found = line[1:].split()[0] == record_id
                continue
            if found:
                seq.append(line.upper())
    dna = "".join(seq).replace(" ", "").replace("\t", "")
    if not dna:
        raise SystemExit(f"ERROR: Transcript not found in FASTA: {record_id}")
    return dna


def translate(dna: str, frame0: int) -> str:
    aa = []
    for i in range(frame0, len(dna) - 2, 3):
        aa.append(CODON_TABLE.get(dna[i:i+3], "X"))
    return "".join(aa)


dna = read_fasta_record(fasta, target)
hits = []
for frame0 in (0, 1, 2):
    aa = translate(dna, frame0)
    idx = aa.find(peptide)
    while idx != -1:
        aa_start = idx + 1
        aa_end = idx + len(peptide)
        nt_start = frame0 + idx * 3 + 1
        nt_end = frame0 + (idx + len(peptide)) * 3
        hits.append((frame0 + 1, aa_start, aa_end, nt_start, nt_end))
        idx = aa.find(peptide, idx + 1)

if not hits:
    raise SystemExit(f"ERROR: Peptide {peptide} not found in {target} in forward-frame translation.")

with out_path.open("w", encoding="utf-8") as out:
    out.write("Transcript\tPeptide\tFrame\tAA_start\tAA_end\tNT_start\tNT_end\tNT_len\n")
    for frame, aa_s, aa_e, nt_s, nt_e in hits:
        out.write(f"{target}\t{peptide}\t{frame}\t{aa_s}\t{aa_e}\t{nt_s}\t{nt_e}\t{nt_e - nt_s + 1}\n")
PY

# 2) Project peptide transcript nucleotide interval to genomic coordinates using CIGAR.
python3 - "${OUTDIR}/EVA_hg38.bam" "${EVA_FASTA}" "${TARGET_TX}" "${PEP}" "${OUTDIR}/peptide_${PEP}_genomic_coords.tsv" <<'PY'
import csv
import re
import subprocess
import sys
from pathlib import Path

bam = Path(sys.argv[1])
fasta = Path(sys.argv[2])
qname = sys.argv[3]
peptide = sys.argv[4]
out_path = Path(sys.argv[5])

hit_file = out_path.parent / f"peptide_in_{qname}_transcript.tsv"
with hit_file.open("r", encoding="utf-8") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    first_hit = next(reader, None)
if first_hit is None:
    raise SystemExit(f"ERROR: No peptide hit rows found in {hit_file}")

nt_start = int(first_hit["NT_start"])
nt_end = int(first_hit["NT_end"])

seq = []
found = False
with fasta.open("r", encoding="utf-8") as handle:
    for line in handle:
        line = line.strip()
        if line.startswith(">"):
            found = line[1:].split()[0] == qname
            continue
        if found:
            seq.append(line.upper())
transcript_len = len("".join(seq))

sam_lines = subprocess.run(["samtools", "view", str(bam)], check=True, text=True, capture_output=True).stdout.splitlines()
matching = [line for line in sam_lines if line.startswith(qname + "\t")]
if not matching:
    raise SystemExit(f"ERROR: No alignment found for {qname} in {bam}")

fields = matching[0].split("\t")
flag = int(fields[1])
chrom = fields[2]
pos = int(fields[3])
cigar = fields[5]
is_reverse = (flag & 16) != 0
strand = "-" if is_reverse else "+"

# Convert transcript coordinates to aligned query coordinates if the transcript
# aligned to the reverse strand.
if is_reverse:
    query_start = transcript_len - nt_end + 1
    query_end = transcript_len - nt_start + 1
else:
    query_start, query_end = nt_start, nt_end

ops = [(int(n), op) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)]
ref_pos = pos
query_pos = 1
intervals = []

for length, op in ops:
    if op in {"M", "=", "X"}:
        q1, q2 = query_pos, query_pos + length - 1
        r1 = ref_pos
        ov1 = max(q1, query_start)
        ov2 = min(q2, query_end)
        if ov1 <= ov2:
            intervals.append((chrom, r1 + (ov1 - q1), r1 + (ov2 - q1), strand))
        query_pos += length
        ref_pos += length
    elif op in {"I", "S"}:
        query_pos += length
    elif op in {"D", "N"}:
        ref_pos += length
    elif op in {"H", "P"}:
        continue

with out_path.open("w", encoding="utf-8") as out:
    out.write("Peptide\tTranscript\tTranscript_NT_start\tTranscript_NT_end\tChrom\tGenomic_start_1based\tGenomic_end_1based\tStrand\tLength_nt\n")
    for chrom_i, start, end, strand_i in intervals:
        out.write(f"{peptide}\t{qname}\t{nt_start}\t{nt_end}\t{chrom_i}\t{start}\t{end}\t{strand_i}\t{end - start + 1}\n")
PY

# 3) Build target transcript locus BED from Script 01 locus summary.
python3 - "${OUTDIR}/EVA_transcript_loci.tsv" "${TARGET_TX}" "${OUTDIR}/${TARGET_TX}_locus_hg38.bed" <<'PY'
import csv
import sys
from pathlib import Path

loci_path = Path(sys.argv[1])
target = sys.argv[2]
out_path = Path(sys.argv[3])

found = False
with loci_path.open("r", encoding="utf-8") as handle, out_path.open("w", encoding="utf-8") as out:
    reader = csv.DictReader(handle, delimiter="\t")
    for row in reader:
        if row["Transcript"] == target:
            chrom = row["Chrom"]
            start0 = int(row["Start_1based"]) - 1
            end0 = int(row["End_1based"])
            strand = row["Strand"]
            out.write(f"{chrom}\t{start0}\t{end0}\t{target}_locus\t0\t{strand}\n")
            found = True
            break
if not found:
    raise SystemExit(f"ERROR: Target transcript not found in locus table: {target}")
PY

bedtools sort -g "${OUTDIR}/hg38.genome" -i "${OUTDIR}/${TARGET_TX}_locus_hg38.bed" \
  > "${OUTDIR}/${TARGET_TX}_locus_hg38.sorted.bed"

# 4) Repeat composition by class. The RepeatMasker class column is configurable
# as a 1-based column index within the RMSK BED file, not the intersect output.
bedtools intersect -wo \
  -a "${OUTDIR}/${TARGET_TX}_locus_hg38.sorted.bed" \
  -b "${RMSK_BED}" \
  > "${OUTDIR}/${TARGET_TX}_repeatmasker_intersections.tsv"

python3 - "${OUTDIR}/${TARGET_TX}_locus_hg38.sorted.bed" "${OUTDIR}/${TARGET_TX}_repeatmasker_intersections.tsv" "${RMSK_CLASS_COL}" "${OUTDIR}/${TARGET_TX}_repeat_composition_clean.tsv" <<'PY'
import sys
from collections import defaultdict
from pathlib import Path

locus_bed = Path(sys.argv[1])
intersections = Path(sys.argv[2])
rmsk_class_col = int(sys.argv[3])  # 1-based, within B/RMSK BED columns
out_path = Path(sys.argv[4])

with locus_bed.open("r", encoding="utf-8") as handle:
    rows = [line.rstrip("\n").split("\t") for line in handle if line.strip()]
if len(rows) != 1:
    raise SystemExit("ERROR: Expected a single target locus interval.")

locus = rows[0]
locus_chrom = locus[0]
locus_start = int(locus[1])
locus_end = int(locus[2])
locus_len = locus_end - locus_start

a_cols = 6
class_to_intervals = defaultdict(list)

if intersections.stat().st_size > 0:
    with intersections.open("r", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= a_cols + rmsk_class_col - 1:
                raise SystemExit("ERROR: RepeatMasker class column index is outside the available columns. Check parameters.rmsk_class_column.")
            b_chrom = fields[a_cols]
            b_start = int(fields[a_cols + 1])
            b_end = int(fields[a_cols + 2])
            repeat_class = fields[a_cols + rmsk_class_col - 1]
            ov_start = max(locus_start, b_start)
            ov_end = min(locus_end, b_end)
            if b_chrom == locus_chrom and ov_end > ov_start:
                class_to_intervals[repeat_class].append((ov_start, ov_end))


def merge_len(intervals):
    if not intervals:
        return 0
    intervals = sorted(intervals)
    merged = []
    cur_s, cur_e = intervals[0]
    for s, e in intervals[1:]:
        if s <= cur_e:
            cur_e = max(cur_e, e)
        else:
            merged.append((cur_s, cur_e))
            cur_s, cur_e = s, e
    merged.append((cur_s, cur_e))
    return sum(e - s for s, e in merged)

covered = {cls: merge_len(intervals) for cls, intervals in class_to_intervals.items()}
annotated_bp = sum(covered.values())
unannotated_bp = max(0, locus_len - annotated_bp)

rows_out = [(cls, bp, 100 * bp / locus_len) for cls, bp in covered.items() if bp > 0]
rows_out.append(("Unannotated", unannotated_bp, 100 * unannotated_bp / locus_len))
rows_out.sort(key=lambda x: x[1], reverse=True)

with out_path.open("w", encoding="utf-8") as out:
    out.write("Repeat_class\tbp_covered\tpercent_of_locus\n")
    for cls, bp, pct in rows_out:
        out.write(f"{cls}\t{bp}\t{pct:.2f}\n")
PY

# 5) LTR-only RepeatMasker hits for reporting. Uses the same configurable class column.
python3 - "${OUTDIR}/${TARGET_TX}_repeatmasker_intersections.tsv" "${RMSK_CLASS_COL}" "${OUTDIR}/${TARGET_TX}_locus_LTR_hits.tsv" <<'PY'
import sys
from pathlib import Path

intersections = Path(sys.argv[1])
rmsk_class_col = int(sys.argv[2])
out_path = Path(sys.argv[3])
a_cols = 6

with out_path.open("w", encoding="utf-8") as out:
    if intersections.stat().st_size > 0:
        with intersections.open("r", encoding="utf-8") as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                repeat_class = fields[a_cols + rmsk_class_col - 1]
                if repeat_class == "LTR":
                    out.write(line)
PY

# 6) Extract the peptide-coding 27-nt sequence from the transcript.
python3 - "${EVA_FASTA}" "${TARGET_TX}" "${PEP}" "${OUTDIR}/${PEP}_27nt.fa" <<'PY'
import csv
import sys
from pathlib import Path

fasta = Path(sys.argv[1])
target = sys.argv[2]
peptide = sys.argv[3]
out_path = Path(sys.argv[4])
hit_path = out_path.parent / f"peptide_in_{target}_transcript.tsv"

with hit_path.open("r", encoding="utf-8") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    first_hit = next(reader, None)
if first_hit is None:
    raise SystemExit(f"ERROR: No peptide hit rows found in {hit_path}")

nt_start = int(first_hit["NT_start"])
nt_end = int(first_hit["NT_end"])

seq = []
found = False
with fasta.open("r", encoding="utf-8") as handle:
    for line in handle:
        line = line.strip()
        if line.startswith(">"):
            found = line[1:].split()[0] == target
            continue
        if found:
            seq.append(line.upper())

dna = "".join(seq)
query = dna[nt_start - 1:nt_end]
if len(query) != len(peptide) * 3:
    raise SystemExit(f"ERROR: Extracted nucleotide length is {len(query)}; expected {len(peptide) * 3}.")

with out_path.open("w", encoding="utf-8") as out:
    out.write(f">{peptide}_{len(query)}nt_transcript\n")
    out.write(query + "\n")
PY

echo "Completed target peptide and repeat-context analysis. Outputs written to: ${OUTDIR}"
