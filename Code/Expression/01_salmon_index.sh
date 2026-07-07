#!/usr/bin/env bash
set -euo pipefail

# Build Salmon index for EVA exon-only reference.
#
# Input:
#   reference/EVA_exons_only.fa
#
# Output:
#   salmon_index/EVA_exons_only_index/
#
# Usage:
#   bash code/expression/01_make_salmon_index.sh
#
# Requirements:
#   salmon

REF_FASTA="${REF_FASTA:-reference/EVA_exons_only.fa}"
INDEX_DIR="${INDEX_DIR:-salmon_index/EVA_exons_only_index}"
KMER_SIZE="${KMER_SIZE:-31}"

if ! command -v salmon >/dev/null 2>&1; then
    echo "ERROR: salmon is not available in PATH." >&2
    exit 1
fi

if [[ ! -f "${REF_FASTA}" ]]; then
    echo "ERROR: Reference FASTA not found: ${REF_FASTA}" >&2
    exit 1
fi

mkdir -p "$(dirname "${INDEX_DIR}")"

echo "Building Salmon index"
echo "Reference: ${REF_FASTA}"
echo "Index:     ${INDEX_DIR}"
echo "K-mer:     ${KMER_SIZE}"

salmon index \
    -t "${REF_FASTA}" \
    -i "${INDEX_DIR}" \
    -k "${KMER_SIZE}"

echo "Done: ${INDEX_DIR}"