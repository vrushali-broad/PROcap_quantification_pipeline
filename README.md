# PRO-Cap Quantification Pipeline

## Author
**Vrushali D. Fangal**  
**Date:** 2024-05-23  

## Overview
The **PRO-cap Quantification Pipeline** is designed for the preprocessing of PRO-cap (or similar) paired-end sequencing data, converting raw FASTQ files into BigWig files for downstream analysis. This pipeline integrates multiple tools to perform quality control, trimming, alignment, and various post-processing steps to generate high-quality transcription start site (TSS) and pause site data.

## Features
- **Quality control**: Runs FastQC on raw reads.
- **Adapter trimming & UMI extraction**: Uses Fastp for trimming and extracting UMIs.
- **Genome handling**: Downloads genome files and annotations, creates STAR indices.
- **Alignment**: Aligns reads to a reference genome using STAR.
- **BAM processing**: Filtering, deduplication, indexing, and replicate merging.
- **BigWig generation**: Creates BedGraph and BigWig files for 5' and 3' strand ends.
- **Optional**: BigWig generation using PINTS Visualizer.

## Workflow
1. **Input Validation & Preprocessing**: Process paired-end FASTQ input files.
2. **Quality Control**: Assess raw read quality.
3. **Trimming & UMI Extraction**: Remove adapters, extract UMIs, and filter reads.
4. **Genome Preparation**:
   - Fetch human genome (hg38) with optional spike-in references.
   - Create STAR genome indices for alignment.
5. **Read Alignment**:
   - Align reads using STAR with customizable parameters.
   - Filter and remove PCR duplicates.
   - Generate strand-specific BAM files for TSS and pause sites.
6. **Replicate Handling**: Merge BAM files if multiple replicates are provided.
7. **BigWig Generation**:
   - Convert BAM to BigWig files for transcription visualization.
   - (Optional) Run PINTS Visualizer for enhanced data exploration.
8. **Cleanup**: Remove intermediate files to save storage.

## Dependencies
- **Python 3.7+**
- **Required Tools**:
  - FastQC, Fastp, STAR, samtools, bedtools, umi_tools, bedGraphToBigWig
- **Python Libraries**:
  - argparse, subprocess, os, logging, requests, tqdm, gzip, shutil, pyfaidx, datetime

### Setting up the environment
To install the dependencies, create and activate the environment:
```sh
conda env create -f environment.yml
conda activate procap_env
```

## Usage
Run the script as follows:
```sh
python PROcap_processing_pipeline.py \
    -f1 read1.fastq.gz -f2 read2.fastq.gz \
    -o output_dir -p sample_prefix -g genome_dir \
    --umi_length 6 --threads 48 --adapter_seq1 TGGAATTCTCGGGTGCCAAGG \
    --adapter_seq2 GATCGTCGGACTGTAGAACTCTGAAC --run_pints_preprocess
```

## Arguments
| Argument | Description |
|----------|-------------|
| `-f1`, `--fastq1` | Input FASTQ file for read 1. |
| `-f2`, `--fastq2` | Input FASTQ file for read 2. |
| `-o`, `--output_dir` | Directory for output files. (Default: `out_dir`) |
| `-g`, `--genome_dir` | Directory for genome reference files. (Default: `genome_dir`) |
| `-p`, `--prefix` | Prefix for output file names. |
| `-u`, `--umi_length` | Length of UMI to extract. (Default: `6`) |
| `-t`, `--threads` | Number of threads to use. (Default: `48`) |
| `--spike_ins` | Include spike-in references (e.g., Drosophila genome). |
| `--run_pints_preprocess` | Run PINTS visualizer on processed BAM files. |

## Outputs
- **BigWig files**: 5' and 3' strand-specific transcription site data.
- **Quality reports**: FastQC reports for raw and trimmed reads.
- **Processed BAM files**: Aligned, filtered, deduplicated, and indexed BAM files.
- **Merged BAM files**: If multiple replicates are provided.
- **STAR indices**: Genome reference indices for alignment.

---

# Genome-wide Correlation Analysis of BigWig Files

## Overview
This script calculates and visualizes genome-wide correlation between two BigWig files. It bins the genomic data, extracts values, and computes a correlation scatter plot. Multi-processing is used to speed up the per-chromosome binning process.

## Features
- Reads chromosome sizes from an input file.
- Extracts signal values from BigWig files in non-overlapping bins.
- Computes Pearson correlation between log-transformed values.
- Generates a scatter plot with a 45-degree reference line.
- Utilizes multiprocessing to speed up processing.

## Usage
Run the script as follows:
```sh
python procap_correlation.py -f1 file1.bw -f2 file2.bw --bin_size 1000 --chrom_sizes_file chrom.sizes
```

## Dependencies
- **Python 3.7+**
- **Required Libraries**:
  - pyBigWig, numpy, matplotlib, argparse, multiprocessing, tqdm, scipy, logging

## Outputs
- **Scatter plot**: Displays correlation between two BigWig files.
- **Pearson correlation coefficient**: Measures the strength of association between signals.


## Notes
- Ensure all required tools and Python libraries are installed before running the script.
- Intermediate files can be removed using the cleanup functionality.

For further details, please contact **Vrushali D. Fangal at YuLab, Cornell University**.

