I'll help you create a comprehensive wrapper script for the LOH-HLA analysis. Based on the context and requirements, here's a detailed implementation:

```R
#!/usr/bin/env Rscript

library(optparse)
library(parallel)

# Define option list matching LOHHLAscript.R parameters
option_list <- list(
  make_option(c("-b", "--bamList"), type="character", default=NULL,
              help="Text file containing paths to BAM files (REQUIRED)"),
  make_option(c("-h", "--hlaList"), type="character", default=NULL,
              help="Text file containing paths to HLA type files (REQUIRED)"),
  make_option(c("-o", "--outputDir"), type="character", default=getwd(),
              help="Output directory for analysis results [default: current working directory]"),
  make_option(c("-n", "--normalBAM"), type="character", default="FALSE",
              help="Normal BAM file path. Can be FALSE to run without normal sample [default: FALSE]"),
  make_option(c("-hl", "--hlaFastaLoc"), type="character", default="~/lohhla/data/hla_all.fasta",
              help="Location of HLA FASTA file [default: %default]"),
  make_option(c("-cn", "--copyNumLoc"), type="character", default="FALSE", 
              help="Location of patient purity and ploidy output. Can be FALSE [default: FALSE]"),
  make_option(c("-mc", "--minCoverage"), type="numeric", default=30,
              help="Minimum coverage at mismatch site [default: %default]"),
  make_option(c("-k", "--kmerSize"), type="numeric", default=50,
              help="Size of kmers to fish with [default: %default]"),
  make_option(c("-mm", "--numMisMatch"), type="numeric", default=1,
              help="Number of mismatches allowed in read to map to HLA allele [default: %default]"),
  make_option(c("--mappingStep"), type="logical", default=TRUE,
              help="Perform mapping to HLA alleles [default: %default]"),
  make_option(c("--fishingStep"), type="logical", default=TRUE,
              help="Look for fished reads matching kmers [default: %default]"),
  make_option(c("--plottingStep"), type="logical", default=TRUE,
              help="Generate plots [default: %default]"),
  make_option(c("--coverageStep"), type="logical", default=TRUE,
              help="Analyze coverage differences [default: %default]"),
  make_option(c("--cleanUp"), type="logical", default=TRUE,
              help="Remove temporary files [default: %default]"),
  make_option(c("--novoDir"), type="character", default="",
              help="Path to novoalign executable [default: empty]"),
  make_option(c("--gatkDir"), type="character", default="",
              help="Path to GATK executable [default: empty]"),
  make_option(c("--hlaExonLoc"), type="character", default="~/lohhla/data/hla.dat",
              help="HLA exon boundaries for plotting [default: %default]"),
  make_option(c("--ignoreWarnings"), type="logical", default=TRUE,
              help="Continue running with warnings [default: %default]"),
  make_option(c("-c", "--cores"), type="numeric", default=1,
              help="Number of cores to use for parallel processing [default: %default]")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate required inputs
if(is.null(opt$bamList) || is.null(opt$hlaList)) {
  print_help(opt_parser)
  stop("BAM list and HLA list are required arguments.", call.=FALSE)
}

# Read input files
bam_files <- readLines(opt$bamList)
hla_files <- readLines(opt$hlaList)

# Validate input lists match
if(length(bam_files) != length(hla_files)) {
  stop("Number of BAM files must match number of HLA files.")
}

# Function to run LOHHLAscript for each sample
run_lohhla <- function(bam_file, hla_file, opt) {
  # Construct command dynamically
  cmd_args <- c(
    "--patientId", basename(bam_file),
    "--outputDir", opt$outputDir,
    "--normalBAMfile", opt$normalBAM,
    "--tumorBAMfile", bam_file,
    "--hlaPath", hla_file,
    "--HLAfastaLoc", opt$hlaFastaLoc,
    "--CopyNumLoc", opt$copyNumLoc,
    "--minCoverageFilter", opt$minCoverage,
    "--kmerSize", opt$kmerSize,
    "--numMisMatch", opt$numMisMatch,
    "--mappingStep", opt$mappingStep,
    "--fishingStep", opt$fishingStep,
    "--plottingStep", opt$plottingStep,
    "--coverageStep", opt$coverageStep,
    "--cleanUp", opt$cleanUp,
    "--novoDir", opt$novoDir,
    "--gatkDir", opt$gatkDir,
    "--HLAexonLoc", opt$hlaExonLoc,
    "--ignoreWarnings", opt$ignoreWarnings
  )
  
  # Run LOHHLAscript
  system2("Rscript", args = c("LOHHLAscript.R", cmd_args))
}

# Parallel execution
mclapply(seq_along(bam_files), function(i) {
  run_lohhla(bam_files[i], hla_files[i], opt)
}, mc.cores = opt$cores)
```

And here's a corresponding README.md:

```markdown
# LOH-HLA Wrapper Script

## Overview
This wrapper script runs Loss of Heterozygosity (LOH) analysis on HLA regions for multiple samples.

## Usage
```bash
Rscript lohhla_wrapper.R -b bam_list.txt -h hla_list.txt [options]
```

## Required Arguments
- `-b, --bamList`: Text file with paths to tumor BAM files
- `-h, --hlaList`: Text file with paths to corresponding HLA type files

## Optional Arguments
### Input/Output
- `-o, --outputDir`: Output directory [default: current working directory]
- `-n, --normalBAM`: Normal BAM file path [default: FALSE]

### HLA Analysis Parameters
- `-hl, --hlaFastaLoc`: HLA FASTA file location
- `-mc, --minCoverage`: Minimum coverage at mismatch site [default: 30]
- `-k, --kmerSize`: K-mer size for fishing [default: 50]
- `-mm, --numMisMatch`: Allowed mismatches in read mapping [default: 1]

### Analysis Steps
- `--mappingStep`: Perform mapping to HLA alleles [default: TRUE]
- `--fishingStep`: Look for fished reads [default: TRUE]
- `--plottingStep`: Generate plots [default: TRUE]
- `--coverageStep`: Analyze coverage differences [default: TRUE]

### Utility Options
- `--cleanUp`: Remove temporary files [default: TRUE]
- `--ignoreWarnings`: Continue running with warnings [default: TRUE]
- `-c, --cores`: Number of parallel cores [default: 1]

### Tool Paths
- `--novoDir`: Path to novoalign executable
- `--gatkDir`: Path to GATK executable
- `--hlaExonLoc`: HLA exon boundaries file
```

Example command line:
```bash
Rscript lohhla_wrapper.R -b tumor_bams.txt -h hla_types.txt -o /path/to/output -c 4
```

This wrapper script:
1. Matches all parameters from the original LOHHLAscript.R
2. Allows parallel processing of multiple samples
3. Validates input files
4. Dynamically constructs command for each sample
5. Provides comprehensive help and documentation

Note: Ensure LOHHLAscript.R is in the same directory or provide full path.