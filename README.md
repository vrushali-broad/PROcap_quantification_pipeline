# PROcap_quantification_pipeline

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

## Outputs
- **BigWig files**: 5' and 3' strand-specific transcription site data.
- **Quality reports**: FastQC reports for raw and trimmed reads.
- **Processed BAM files**: Aligned, filtered, deduplicated, and indexed BAM files.
- **Merged BAM files**: If multiple replicates are provided.
- **STAR indices**: Genome reference indices for alignment.

## Notes
- Ensure all required tools and Python libraries are installed before running the script.
- Intermediate files can be removed using the cleanup functionality.

For further details, please contact **YuLab, Cornell University**.

