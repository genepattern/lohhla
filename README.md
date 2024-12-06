# LOH-HLA Analysis Wrapper Script

## Description
A comprehensive bioinformatics tool for Loss of Heterozygosity (LOH) analysis in HLA regions, designed for parallel processing of multiple genomic samples.

## Module Details
- **Authors**: Edwin Huang + LLMs
- **Categories**: HLA analysis
- **Source Repository**: [https://github.com/mskcc/lohhla](https://github.com/mskcc/lohhla)
- **Contact**: edh021@cloud.ucsd.edu

## Input Files
1. **BAM List File** (Text file)
   - Contains paths to tumor BAM files
2. **HLA Type List File** (Text file)
   - Contains corresponding HLA type file paths
3. **Optional: Normal BAM file**

## Output Files
1. Coverage plots
2. HLA allele mapping results
3. LOH analysis reports
4. Temporary intermediate files (if clean-up disabled)

## Parameters

| Parameter | Description | Default Value | Type |
|-----------|-------------|---------------|------|
| Output Directory | Specifies location for analysis results | Current working directory | File path |
| Normal BAM File | Optional reference normal sample | FALSE | File path |
| HLA FASTA Location | Reference HLA sequence database | ~/lohhla/data/hla_all.fasta | File path |
| HLA Exon Location | HLA exon boundary information for plotting | ~/lohhla/data/hla.dat | File path |
| Minimum Coverage | Minimum read coverage at mismatch sites | 30 | Numeric |
| K-mer Size | Size of genomic fragments for read mapping | 50 | Numeric |
| Mismatch Tolerance | Maximum allowed mismatches in read-to-allele mapping | 1 | Numeric |
| Mapping Step | Perform mapping to HLA alleles | TRUE | Boolean |
| Fishing Step | Identify reads matching specific k-mers | TRUE | Boolean |
| Plotting Step | Generate visualization of analysis results | TRUE | Boolean |
| Coverage Step | Analyze coverage differences across regions | TRUE | Boolean |
| Clean Up | Remove temporary files after analysis | TRUE | Boolean |
| Ignore Warnings | Continue execution despite non-critical warnings | TRUE | Boolean |
| Parallel Cores | Number of CPU cores for parallel processing | 1 | Numeric |
| Novoalign Directory | Path to alignment tool executable | Empty | File path |
| GATK Directory | Path to Genome Analysis Toolkit executable | Empty | File path |

## Computational Requirements
- R environment
- Required R libraries: optparse, parallel
- External tools: Novoalign, GATK (optional)
- Sufficient computational resources based on sample complexity and selected cores

## Use Cases
- Cancer genomics research
- HLA region variation analysis
- Immunogenomics studies
- Personalized medicine investigations

## Limitations
- Requires matched BAM and HLA type files
- Performance dependent on input data quality
- Computational intensity increases with sample complexity

## Recommended Workflow
1. Prepare input BAM and HLA type files
2. Configure analysis parameters
3. Select appropriate computational resources
4. Execute wrapper script
5. Review generated results