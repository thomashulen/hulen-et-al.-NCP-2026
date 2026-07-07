#!/usr/bin/env bash
set -euo pipefail

# Quantify EVA exon-level expression using Salmon.
#
# Inputs:
#   metadata/sample_sheet.tsv
#   salmon_index/EVA_exons_only_index/
#   FASTQ files listed in the sample sheet
#
# Output:
#   results/salmon_exon_quant/<sample>/quant.sf
#
# Usage:
#   bash code/expression/02_quantify_exons.sh
#
# Requirements:
#   salmon
#
# Expected sample sheet columns:
#   sample    group    fastq_1    fastq_2

SAMPLE_SHEET="${SAMPLE_SHEET:-metadata/sample_sheet.tsv}"
INDEX_DIR="${INDEX_DIR:-salmon_index/EVA_exons_only_index}"
OUT_DIR="${OUT_DIR:-results/salmon_exon_quant}"

LIB_TYPE="${LIB_TYPE:-A}"
THREADS="${THREADS:-8}"

if ! command -v salmon >/dev/null 2>&1; then
    echo "ERROR: salmon is not available in PATH." >&2
    exit 1
fi

if [[ ! -f "${SAMPLE_SHEET}" ]]; then
    echo "ERROR: Sample sheet not found: ${SAMPLE_SHEET}" >&2
    exit 1
fi

if [[ ! -d "${INDEX_DIR}" ]]; then
    echo "ERROR: Salmon index not found: ${INDEX_DIR}" >&2
    echo "Run code/expression/01_make_salmon_index.sh first." >&2
    exit 1
fi

mkdir -p "${OUT_DIR}"

echo "Starting Salmon quantification"
echo "Sample sheet: ${SAMPLE_SHEET}"
echo "Index:        ${INDEX_DIR}"
echo "Output dir:   ${OUT_DIR}"

tail -n +2 "${SAMPLE_SHEET}" | while IFS=$'\t' read -r sample group fastq_1 fastq_2; do
    [[ -z "${sample}" ]] && continue

    echo "------------------------------------------------------------"
    echo "Processing sample: ${sample}"
    echo "Group:             ${group}"

    if [[ ! -f "${fastq_1}" ]]; then
        echo "ERROR: FASTQ 1 not found for sample ${sample}: ${fastq_1}" >&2
        exit 1
    fi

    if [[ ! -f "${fastq_2}" ]]; then
        echo "ERROR: FASTQ 2 not found for sample ${sample}: ${fastq_2}" >&2
        exit 1
    fi

    salmon quant \
        -i "${INDEX_DIR}" \
        -l "${LIB_TYPE}" \
        -1 "${fastq_1}" \
        -2 "${fastq_2}" \
        -p "${THREADS}" \
        --validateMappings \
        -o "${OUT_DIR}/${sample}"

    echo "Finished sample: ${sample}"
done

echo "Done. Salmon quantification results written to: ${OUT_DIR}"