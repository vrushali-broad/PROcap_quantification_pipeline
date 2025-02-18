# PROcap_quantification_pipeline

**Author:** Vrushali D. Fangal  
**Date:** 2024-05-23  
**Organization:** YuLab, Cornell University  

## Description  
This script provides a comprehensive preprocessing pipeline for PRO-cap (or similar) paired-end sequencing data, converting raw FASTQ files into BigWig files for downstream analysis. It integrates multiple tools and performs quality control, trimming, alignment, and various post-processing steps.

## Key Features  
- Quality control with **FastQC**  
- Adapter trimming and UMI extraction using **Fastp**  
- Genome download, indexing, and annotation retrieval  
- Alignment of reads to a reference genome using **STAR**  
- BAM file processing: filtering, deduplication, and indexing  
- Handling of multiple replicates with automatic merging of BAM files  
- Generation of **BedGraph** and **BigWig** files for **5' and 3'** strand ends  
- Optional **BigWig generation using PINTS Visualizer**  

## Pipeline Overview  
1. **Validate and preprocess** paired-end FASTQ input files.  
2. **Perform quality control** to assess raw read quality.  
3. **Trim adapters, extract UMIs, and filter reads.**  
4. **Download genome and annotations** for alignment:  
   - Fetch human genome (**hg38**) by default, with optional spike-in references.  
   - Create **STAR** genome indices for alignment.  
5. **Align reads using STAR** with customizable parameters.  
6. **Process alignment outputs:**  
   - Filter chromosomes of interest.  
   - Remove PCR duplicates using **UMI-tools**.  
   - Separate and generate **strand-specific BAM files** for **read1 (3' pause site) and read2 (5' initiation site)**.  
7. **Merge replicates** into a single BAM file for downstream analysis.  
8. **Generate BedGraph and BigWig files** for visualization of transcription start sites (**TSS**) and pause sites.  
9. **(Optional) Run PINTS Visualizer** for additional data exploration.  
10. **Cleanup intermediate files** to save storage.  

## Dependencies  
- **Python 3.7+**  
- **External tools:**  
  - FastQC  
  - Fastp  
  - STAR  
  - samtools  
  - bedtools  
  - umi_tools  
  - bedGraphToBigWig  
- **Python libraries:**  
  - `argparse`  
  - `subprocess`  
  - `os`  
  - `logging`  
  - `requests`  
  - `tqdm`  
  - `gzip`  
  - `shutil`  
  - `pyfaidx`  
  - `datetime`  


## Installation

### Setting Up the Conda Environment
To set up the environment using Conda, follow these steps:

```bash
conda env create -f environment.yml
```

This will create a new Conda environment with all the necessary dependencies.

### Activating the Environment
Once the installation is complete, activate the environment:

```bash
conda activate procap_env
```

Ensure that `procap_env` is the name specified in your `environment.yml` file.

If you need to update the environment after modifying the `environment.yml` file, run:

```bash
conda env update --file environment.yml --prune
```

### Installing Python Dependencies (Alternative)
If using Conda is not an option, you can install the required Python dependencies manually using `pip`:

```bash
pip install tqdm pyfaidx requests
```

## Usage

Run the pipeline with the following command:

```bash
python PROcap_processing_pipeline.py \
    -f1 read1.fastq.gz -f2 read2.fastq.gz \
    -o output_dir -p sample_prefix -g genome_dir \
    --umi_length 6 --threads 48 --adapter_seq1 TGGAATTCTCGGGTGCCAAGG \
    --adapter_seq2 GATCGTCGGACTGTAGAACTCTGAAC --run_pints
```

## Arguments

| Argument | Description |
|----------|-------------|
| `-f1, --fastq1` | Input FASTQ file for read 1 |
| `-f2, --fastq2` | Input FASTQ file for read 2 |
| `-o, --output_dir` | Directory for output files (default: `out_dir`) |
| `-g, --genome_dir` | Directory for genome reference files (default: `genome_dir`) |
| `-p, --prefix` | Prefix for output file names |
| `-u, --umi_length` | Length of UMI to extract (default: `6`) |
| `-t, --threads` | Number of threads to use (default: `48`) |
| `--spike_ins` | Include spike-in references (e.g., Drosophila genome) |
| `--run_pints_preprocess` | Run PINTS visualizer on processed BAM files |

## Outputs
- BigWig files for 5' and 3' strand ends
- Quality reports, trimmed FASTQ files, and aligned BAM files
- Merged BAM files for replicate handling
- STAR indices, deduplicated BAM files, and filtered outputs

## Example Workflow

```bash
# Run pipeline with 2 paired-end FASTQ files
python PROcap_processing_pipeline.py -f1 sample_1.fastq.gz sample_2.fastq.gz \
    -f2 sample_1.fastq.gz sample_2.fastq.gz -p experiment1 \
    --genome_dir genome/ --threads 16
```

## Notes
- Ensure all required tools and Python libraries are installed before running the script.
- Intermediate files can be removed using the cleanup functionality.



## Installation  
Make sure you have all dependencies installed before running the pipeline.  

Using **conda**, install external dependencies:  
```bash
conda install -c bioconda fastqc fastp star samtools bedtools umi_tools
```

## Install Python Dependencies Using Pip

```bash
pip install tqdm pyfaidx requests
```


# PRO-cap Quantification Pipeline

## Installation

Install Python dependencies using pip:

```bash
pip install tqdm pyfaidx requests
```

## Usage

Run the pipeline with the following command:

```bash
python PROcap_processing_pipeline.py \
    -f1 read1.fastq.gz -f2 read2.fastq.gz \
    -o output_dir -p sample_prefix -g genome_dir \
    --umi_length 6 --threads 48 --adapter_seq1 TGGAATTCTCGGGTGCCAAGG \
    --adapter_seq2 GATCGTCGGACTGTAGAACTCTGAAC --run_pints
```

## Arguments

| Argument | Description |
|----------|-------------|
| `-f1, --fastq1` | Input FASTQ file for read 1 |
| `-f2, --fastq2` | Input FASTQ file for read 2 |
| `-o, --output_dir` | Directory for output files (default: `out_dir`) |
| `-g, --genome_dir` | Directory for genome reference files (default: `genome_dir`) |
| `-p, --prefix` | Prefix for output file names |
| `-u, --umi_length` | Length of UMI to extract (default: `6`) |
| `-t, --threads` | Number of threads to use (default: `48`) |
| `--spike_ins` | Include spike-in references (e.g., Drosophila genome) |
| `--run_pints_preprocess` | Run PINTS visualizer on processed BAM files |

## Outputs

- BigWig files for 5' and 3' strand ends
- Quality reports, trimmed FASTQ files, and aligned BAM files
- Merged BAM files for replicate handling
- STAR indices, deduplicated BAM files, and filtered outputs

## Example Workflow

```bash
# Run pipeline with 2 paired-end fastq files
python PROcap_processing_pipeline.py -f1 sample_1.fastq.gz sample_2.fastq.gz \
    -f2 sample_1.fastq.gz sample_2.fastq.gz -p experiment1 \
    --genome_dir genome/ --threads 16
```

## Notes

- Ensure all required tools and Python libraries are installed before running the script.
- Intermediate files can be removed using the cleanup functionality.






# PRO-cap Quantification Pipeline

## Author
**Vrushali D. Fangal**  
**Date:** 2024-05-23  
**Organization:** YuLab, Cornell University  

## Installation

### Setting Up the Conda Environment
To set up the environment using Conda, follow these steps:

```bash
conda env create -f environment.yml
```

This will create a new Conda environment with all the necessary dependencies.

### Activating the Environment
Once the installation is complete, activate the environment:

```bash
conda activate procap_env
```

Ensure that `procap_env` is the name specified in your `environment.yml` file.

If you need to update the environment after modifying the `environment.yml` file, run:

```bash
conda env update --file environment.yml --prune
```

### Installing Python Dependencies (Alternative)
If using Conda is not an option, you can install the required Python dependencies manually using `pip`:

```bash
pip install tqdm pyfaidx requests
```

## Usage

Run the pipeline with the following command:

```bash
python PROcap_processing_pipeline.py \
    -f1 read1.fastq.gz -f2 read2.fastq.gz \
    -o output_dir -p sample_prefix -g genome_dir \
    --umi_length 6 --threads 8 \
    --adapter_seq1 TGGAATTCTCGGGTGCCAAGG \
    --adapter_seq2 GATCGTCGGACTGTAGAACTCTGAAC \
    --run_pints
```

### Running the Genome Correlation Script
To compute and visualize the correlation between two BigWig files, use the following command:

```bash
python procap_correlation.py \
    -f1 bigwig_file1.bw -f2 bigwig_file2.bw \
    --chrom_sizes_file genome.chrom.sizes \
    --bin_size 1000
```

## Arguments

| Argument | Description |
|----------|-------------|
| `-f1, --fastq1` | Input FASTQ file for read 1 |
| `-f2, --fastq2` | Input FASTQ file for read 2 |
| `-o, --output_dir` | Directory for output files (default: `out_dir`) |
| `-g, --genome_dir` | Directory for genome reference files (default: `genome_dir`) |
| `-p, --prefix` | Prefix for output file names |
| `-u, --umi_length` | Length of UMI to extract (default: `6`) |
| `-t, --threads` | Number of threads to use (default: `48`) |
| `--spike_ins` | Include spike-in references (e.g., Drosophila genome) |
| `--run_pints_preprocess` | Run PINTS visualizer on processed BAM files |

## Outputs
- BigWig files for 5' and 3' strand ends
- Quality reports, trimmed FASTQ files, and aligned BAM files
- Merged BAM files for replicate handling
- STAR indices, deduplicated BAM files, and filtered outputs

## Example Workflow

```bash
# Run pipeline with 2 paired-end FASTQ files
python PROcap_processing_pipeline.py -f1 sample_1.fastq.gz sample_2.fastq.gz \
    -f2 sample_1.fastq.gz sample_2.fastq.gz -p experiment1 \
    --genome_dir genome/ --threads 16
```

```bash
# Run genome correlation analysis
python procap_correlation.py \
    -f1 bigwig1.bw -f2 bigwig2.bw \
    --chrom_sizes_file genome.chrom.sizes \
    --bin_size 1000
```

## Notes
- Ensure all required tools and Python libraries are installed before running the script.
- Intermediate files can be removed using the cleanup functionality.

