# EVA003 Genomic Characterization

Pipeline to:
- Map EVA003 transcript to hg38 (minimap2)
- Derive exon structure
- Intersect with RepeatMasker annotations
- Identify nearest annotated genes (e.g., LINC00518)

## Requirements
- minimap2
- samtools
- bedtools
- python3

## Usage
bash code/eva003_genomic_characterization.sh