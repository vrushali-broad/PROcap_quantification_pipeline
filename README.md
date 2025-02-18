# **PRO-Cap Analysis Pipeline**

## **Author**
**Vrushali D. Fangal**  
**Organization:** YuLab, Cornell University  
**Date:** 2024-05-23  

---

## **Overview**
This repository contains two major components:

1. **PRO-Cap Quantification Pipeline** - A pipeline for processing **PRO-cap sequencing data**, performing alignment, filtering, and BigWig generation for visualization of transcription start sites (TSS) and pause sites.
2. **Genome-wide Correlation Analysis** - A script for comparing **two BigWig files**, computing genome-wide correlations, and generating scatter plots.

These scripts integrate multiple bioinformatics tools and are optimized for efficient data processing.

---

## **Table of Contents**
- [Features](#features)
- [Installation](#installation)
- [Pipeline Workflow](#pipeline-workflow)
- [Usage](#usage)
- [Arguments](#arguments)
- [Outputs](#outputs)
- [Example Commands](#example-commands)
- [Genome-wide Correlation Analysis](#genome-wide-correlation-analysis)
- [Environment Setup](#environment-setup)
- [Dependencies](#dependencies)
- [License](#license)
- [Contact](#contact)

---

## **Features**
### **PRO-Cap Quantification Pipeline**
- **Quality control**: Runs **FastQC** on raw reads.
- **Adapter trimming & UMI extraction**: Uses **Fastp** for adapter removal.
- **Genome setup**: Downloads genome files and annotations, builds **STAR indices**.
- **Alignment**: Uses **STAR** for read alignment.
- **BAM processing**: Performs filtering, deduplication, and indexing.
- **BigWig generation**: Produces **strand-specific** BigWig files.
- **PINTS Visualizer** (optional): Creates additional BigWig files for data exploration.

### **Genome-wide Correlation Analysis**
- **Reads chromosome sizes** from an input file.
- **Extracts signal values** in non-overlapping bins from BigWig files.
- **Computes Pearson correlation** on log-transformed values.
- **Generates scatter plots** with a 45-degree reference line.
- **Utilizes multiprocessing** for fast data processing.

---

## **Installation**
### **1. Clone the Repository**
```sh
git clone https://github.com/vrushali-broad/PROcap_quantification_pipeline.git
cd PROcap_quantification_pipeline
```

### **2. Set Up Conda Environment**
To install dependencies:
```sh
conda env create -f environment.yml
conda activate procap_preprocessing
```
To update the environment later:
```sh
conda env export > environment.yml
```

---

## **Pipeline Workflow**
1. **Preprocess paired-end FASTQ files**.
2. **Perform quality control** using **FastQC**.
3. **Trim adapters and extract UMIs**.
4. **Download genome and annotations**, create **STAR indices**.
5. **Align reads to the reference genome** using **STAR**.
6. **Filter and deduplicate BAM files**.
7. **Merge replicates into a single BAM file** (if multiple).
8. **Generate BigWig files for TSS and pause sites**.
9. **(Optional) Run PINTS Visualizer for further analysis**.
10. **Clean up intermediate files**.

---

## **Usage**
Run the pipeline with:
```sh
python PROcap_processing_pipeline.py \
    -f1 read1.fastq.gz -f2 read2.fastq.gz \
    -o output_dir -p sample_prefix -g genome_dir \
    --umi_length 6 --threads 48 --adapter_seq1 TGGAATTCTCGGGTGCCAAGG \
    --adapter_seq2 GATCGTCGGACTGTAGAACTCTGAAC --run_pints_preprocess
```

---

## **Arguments**
| Argument | Description |
|----------|-------------|
| `-f1`, `--fastq1` | Input FASTQ file for read 1. |
| `-f2`, `--fastq2` | Input FASTQ file for read 2. |
| `-o`, `--output_dir` | Directory for output files. (Default: `out_dir`) |
| `-g`, `--genome_dir` | Directory for genome reference files. |
| `-p`, `--prefix` | Prefix for output file names. |
| `-u`, `--umi_length` | Length of UMI to extract (Default: `6`). |
| `-t`, `--threads` | Number of threads to use (Default: `48`). |
| `--spike_ins` | Include spike-in references (e.g., Drosophila genome). |
| `--run_pints_preprocess` | Run PINTS Visualizer for processed BAM files. |

---

## **Outputs**
- **BigWig files**: 5' and 3' strand-specific transcription site data.
- **Quality reports**: FastQC reports for raw and trimmed reads.
- **Processed BAM files**: Aligned, filtered, deduplicated, and indexed BAM files.
- **Merged BAM files**: If multiple replicates are provided.
- **STAR indices**: Genome reference indices for alignment.

---
## **Example Commands**
```sh
python procap_correlation.py -f1 file1.bw -f2 file2.bw --bin_size 1000 --chrom_sizes_file chrom.sizes
```

---

## **Environment Setup**
To recreate the **exact environment**, use the provided `environment.yml`:
```sh
conda env create -f environment.yml
conda activate procap_preprocessing
```

To export the environment after updates:
```sh
conda env export > environment.yml
```

---

## **Dependencies**
- **Python 3.7+**
- **Required Tools**:
  - FastQC, Fastp, STAR, samtools, bedtools, umi_tools, bedGraphToBigWig, pyBigWig
- **Python Libraries**:
  - argparse, subprocess, os, logging, requests, tqdm, gzip, shutil, pyfaidx, numpy, scipy, matplotlib

---

## **License**
This project is licensed under the **MIT License**.

---

## **Contact**
For questions, reach out to **Vrushali D. Fangal** at **YuLab, Cornell University**.

---




