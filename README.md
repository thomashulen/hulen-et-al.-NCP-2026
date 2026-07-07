# Hulen et al. 2026 Computational Analysis

This repository contains the computational analyses and scripts accompanying the manuscript:

> **Receptor-defined targeting of a genomically unique melanoma-enriched noncanonical antigen**

The repository reproduces the genomic characterization, genome-wide uniqueness analyses, and RNA-sequencing quantification performed for the EVA noncanonical peptide (NCP) transcripts.

---

## Repository structure

```text
Code/
    Expression/
        01_make_salmon_index.sh
        02_quantify_exons.sh

    Genomic/
        01_EVA_mapping_to_hg38.sh
        02_EVA003_peptide_repeat_context.sh
        03_27nt_exact_genome_scan.py
        04_scan_27nt_mismatches.py

Data/
    Dataset1_EVA001.fa
    Dataset2_EVA002.fa
    Dataset3_EVA003.fa
    Dataset4_EVA003_600-nt_core.fa
    Dataset5_LLMRPLRIK.fa
    EVA_exons.fa

Metadata/
    sample_sheet.tsv
    Public_Raw_RNA_seq.txt

Results/
    EVA_transcript_loci.tsv
    EVA001_locus_hg38.bed
    EVA002_locus_hg38.bed
    EVA003_locus_hg38.bed
    EVA003_LLMRPLRIK_hg38.bed
    EVA003_repeat_composition_clean.tsv
    EVA_exons_hg38.stranded.clean.bed
```

---

## Repository contents

The repository contains scripts used to:

- map EVA transcripts to the human reference genome (hg38)
- determine exon structure
- identify genomic loci
- characterize repeat composition using RepeatMasker annotations
- identify the genomic position of the LLMRPLRIK peptide-coding sequence
- evaluate genome-wide uniqueness of the 27-nt coding sequence
- quantify exon-level expression from public RNA sequencing datasets using Salmon

---

## Required software

The analyses require:

- Python ≥3.10
- R ≥4.3
- Salmon
- minimap2
- samtools
- bedtools

Python packages:

- PyYAML

---

## Required reference files

The following reference resources are required but are **not included** in this repository:

- hg38 reference genome FASTA
- hg38 FASTA index (.fai)
- GENCODE gene annotation
- RepeatMasker annotation for hg38

The repository contains a template configuration file indicating where these files should be placed.

---

## Analysis workflow

The analyses are intended to be run in the following order.

### 1. Genome mapping

`01_EVA_mapping_to_hg38.sh`

Maps EVA transcripts to the hg38 reference genome and generates transcript loci.

Outputs include:

- `EVA_transcript_loci.tsv`
- `EVA001_locus_hg38.bed`
- `EVA002_locus_hg38.bed`
- `EVA003_locus_hg38.bed`

---

### 2. Genomic context

`02_EVA003_peptide_repeat_context.sh`

Identifies the genomic location of the peptide-coding sequence and characterizes the surrounding genomic context.

Outputs include:
- `EVA003_LLMRPLRIK_hg38.bed`
- `EVA003_repeat_composition_clean.tsv`
---

### 3. Genome-wide uniqueness

`03_27nt_exact_genome_scan.py`

Determines the number of exact genomic matches for the 27-nt peptide-coding sequence.

`04_scan_27nt_mismatches.py`

Performs genome-wide searches allowing a user-defined number of mismatches.

---

### 4. Expression analysis

`01_make_salmon_index.sh`

Builds a Salmon index from the EVA exon reference FASTA.

`02_quantify_exons.sh`

Quantifies exon-level expression from RNA-seq FASTQ files.

---

## Data

Public RNA sequencing datasets analysed in this study are listed in:

```text

Metadata/Public_Raw_RNA_seq.txt

```

---

# Citation

If you use this repository, please cite the accompanying manuscript:

> **Receptor-defined targeting of a genomically unique melanoma-enriched noncanonical antigen**
---
# License

This repository is distributed under the MIT License.
