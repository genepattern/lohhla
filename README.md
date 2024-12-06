# LOHHLA (Loss of Heterozygosity in Human Leukocyte Antigen) Analysis Wrapper

## Description
A bioinformatics tool designed to analyze Loss of Heterozygosity (LOH) in Human Leukocyte Antigen (HLA) regions using sequencing data.

## Module Details
- **Authors**: Edwin Huang + LLMs
- **Categories**: HLA analysis
- **Source Repo**: [LOHHLA GitHub Repository](https://github.com/mskcc/lohhla/blob/master/LOHHLAscript.R)
- **Contact**: edh021@cloud.ucsd.edu

## Input Files
1. **BAM File List**
   - Type: Text file
   - Description: Contains paths to BAM alignment files
   - Required: Yes
   - Format: One BAM file path per line

2. **HLA File List**
   - Type: Text file
   - Description: Contains paths to HLA-related files
   - Required: Yes
   - Format: One HLA file path per line

## Output Files
1. **RDS Results File**: `lohhla_results.rds`
2. **Optional PDF**: `intermediate_results.pdf`

## Parameters

| Parameter | Description | Default Value | Type |
|-----------|-------------|---------------|------|
| `--output_dir (-o)` | Destination for storing analysis results | "./lohhla_output" | String |
| `--plot_intermediate (-p)` | Generate intermediate visualization plots | FALSE | Boolean |
| `--kmer_length` | Length of genomic subsequences for read analysis | 35 | Integer |
| `--min_mapping_quality` | Threshold for read mapping confidence | 20 | Integer |
| `--baf_threshold` | Cutoff for allele frequency variation detection | 0.2 | Decimal |
| `--max_mismatches` | Maximum allowed nucleotide differences in read alignment | 2 | Integer |

## Analysis Workflow
1. Input Validation
2. Output Directory Creation
3. K-mer File Generation
4. Read Mapping Analysis
5. FASTQ File Combination
6. Event Counting
7. Loss of Heterozygosity Calculation
8. Optional Intermediate Result Visualization

## Computational Requirements
- R Environment
- Required R Libraries: optparse
- Recommended System Resources:
  * Multi-core processor
  * Minimum 16GB RAM
  * Sufficient disk space for BAM/FASTQ files

## Typical Use Cases
- HLA region genetic variation analysis
- Cancer immunogenetics research
- Personalized medicine studies

## Execution Example
```bash
Rscript lohhla_wrapper.R \
  -b bam_file_list.txt \
  -h hla_file_list.txt \
  -o ./output_directory \
  -p \
  --kmer_length 40 \
  --min_mapping_quality 25 \
  --baf_threshold 0.3
```

## Limitations
- Requires high-quality sequencing data
- Performance depends on input file quality
- Placeholder functions need domain-specific implementations