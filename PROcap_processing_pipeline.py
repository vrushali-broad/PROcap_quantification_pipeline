"""
PRO-cap Quantification Pipeline

Author: Vrushali D. Fangal  
Date: 2024-05-23  
Organization: YuLab, Cornell University  

Description:  
This script provides a comprehensive preprocessing pipeline for PRO-cap (or similar) paired-end sequencing data, converting raw FASTQ files into BigWig files for downstream analysis. It integrates multiple tools and performs quality control, trimming, alignment, and various post-processing steps.

Key Features:
- Quality control with FastQC.
- Adapter trimming and UMI extraction using Fastp.
- Genome download, indexing, and annotation retrieval.
- Alignment of reads to a reference genome using STAR.
- BAM file processing: filtering, deduplication, and indexing.
- Handling of multiple replicates with automatic merging of BAM files.
- Generation of BedGraph and BigWig files for 5' and 3' strand ends.
- Normalization support for BedGraph files (CPM, RPM, BPM, RPKM).
- RPKM normalization using gene feature lengths from a BED file.
- Optional bigwig generation using PINTS Visualizer.

Pipeline Overview:
1. Validate and preprocess paired-end FASTQ input files.
2. Perform quality control to assess raw read quality.
3. Trim adapters, extract UMIs, and filter reads.
4. Download genome and annotations for alignment:
   - Fetch human genome (hg38) by default, with optional spike-in references.
   - Create STAR genome indices for alignment.
5. Align reads using STAR with customizable parameters.
6. Process alignment outputs:
   - Filter chromosomes of interest.
   - Remove PCR duplicates using UMI-tools.
   - Separate and generate strand-specific BAM files for read1 (3') and read2 (5').
7. Merge replicates into a single BAM file for downstream analysis.
8. Generate BedGraph files and normalize signals:
   - CPM (Counts Per Million), RPM (Reads Per Million), BPM (Bins Per Million).
   - RPKM (Reads Per Kilobase per Million) using feature lengths from a BED file.
9. Convert normalized BedGraph files to BigWig format.
10. (Optional) Run PINTS Visualizer for additional data exploration.
11. Cleanup intermediate files to save storage.

Main Functions:
- setup_logging: Sets up logging to file with timestamp.
- is_gzipped: Checks if a file is gzipped.
- run_command: Executes a shell command and handles errors.
- run_pints_visualizer: Runs the pints_visualizer tool with specified parameters.
- quality_check: Runs FastQC for quality reports on input files.
- fastp_trim_umi_extraction: Trims reads and extracts UMIs using Fastp.
- download_latest_gencode_annotations: Downloads and unzips the latest GENCODE annotations.
- download_genome: Downloads and unzips a genome reference file.
- download_rDNA_curl: Downloads an rDNA sequence from NCBI.
- generate_star_index: Generates a STAR genome index.
- run_star_alignment: Aligns reads with STAR using specified parameters.
- index_bam: Indexes a BAM file using samtools.
- filter_alignments: Filters a BAM file based on chromosomes and quality.
- remove_pcr_duplicates: Removes PCR duplicates using UMI-tools.
- generate_end_bams: Generates BAM files for read1 (pause site) and read2 (initiation) from paired-end reads.
- fetch_chromosome_sizes: Fetches chromosome sizes from a genome FASTA file.
- normalize_bedgraph: Normalizes BedGraph data (CPM, RPM, BPM, RPKM).
- parse_bed_file: Reads feature lengths from a BED file for RPKM normalization.
- safe_file_move: Safely moves files from source to destination.
- generate_bedgraph_and_bigwig: Generates BedGraph and BigWig files from BAM files.
- cleanup_files: Removes specified files from the filesystem.
- count_reads_fastq: Counts the number of reads in a FASTQ file.
- count_reads_bam: Counts the number of reads in a BAM file.
- parse_args: Parses command line arguments.

Dependencies:
- Python 3.7+
- External tools: FastQC, Fastp, STAR, samtools, bedtools, umi_tools, bedGraphToBigWig
- Python libraries: argparse, subprocess, os, logging, requests, tqdm, gzip, shutil, pyfaidx, datetime, numpy

Usage:
    python PROcap_processing_pipeline.py \
        -f1 read1.fastq.gz -f2 read2.fastq.gz \
        -o output_dir -p sample_prefix -g genome_dir \
        --umi_length 6 --threads 48 --adapter_seq1 TGGAATTCTCGGGTGCCAAGG \
        --adapter_seq2 GATCGTCGGACTGTAGAACTCTGAAC --normalization RPKM \
        --features_file genes.bed --run_pints_preprocess

Arguments:
    -f1, --fastq1        Input FASTQ file for read 1.
    -f2, --fastq2        Input FASTQ file for read 2.
    -o, --output_dir     Directory for output files (default: 'out_dir').
    -g, --genome_dir     Directory for genome reference files (default: 'genome_dir').
    -p, --prefix         Prefix for output file names.
    -u, --umi_length     Length of UMI to extract (default: 6).
    -t, --threads        Number of threads to use (default: 48).
    --spike_ins          Include spike-in references (e.g., Drosophila genome).
    --normalization      Normalization method for BedGraph files (CPM, RPM, BPM, RPKM) (default: CPM).
    --features_file      Path to a BED file for RPKM normalization (required for RPKM).
    --run_pints_preprocess         Run PINTS visualizer on processed BAM files.

Outputs:
- BigWig files for 5' and 3' strand ends.
- Normalized BedGraph files (CPM, RPM, BPM, or RPKM).
- Quality reports, trimmed FASTQ files, and aligned BAM files.
- Merged BAM files for replicate handling.
- STAR indices, deduplicated BAM files, and filtered outputs.

Note:
- Ensure all required tools and Python libraries are installed before running the script.
- Intermediate files can be removed using the cleanup functionality.
"""

import argparse
import subprocess
import os
import requests
import logging
import numpy as np
from tqdm import tqdm
import re
import gzip
import shutil
import pyfaidx
from datetime import datetime
import tempfile
import shutil
import sys


################ PRO-cap paired-end reads ##################
######### TSS- 5' ---------------3'                  #######
######### ------->                     Read 2        #######
######### =========================    DNA fragment  #######
#########                <---------    Read 1        #######
############################################################

def setup_logging(output_dir):
    """
    Set up logging to file in the specified output directory with a timestamp.

    Parameters:
    output_dir (str): Path to the directory where the log file will be stored.
    """
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Define the format of the log messages with "EnhancerNet" label before the timestamp
    log_format = 'EnhancerNet - %(asctime)s - %(levelname)s - %(message)s'
    
    # Get the current date and time to append to the log filename
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    log_filename = f"log_{timestamp}.log"

    # Create a full path for the log file
    log_file_path = os.path.join(output_dir, log_filename)

    # Create handlers
    file_handler = logging.FileHandler(log_file_path, mode='w')
    console_handler = logging.StreamHandler()

    # Create formatter and add it to handlers
    formatter = logging.Formatter(log_format, datefmt='%Y-%m-%d %H:%M:%S')
    file_handler.setFormatter(formatter)
    console_handler.setFormatter(formatter)

    # Get the root logger and set level
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    # Add handlers to the logger
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    # Inform the user where log messages will be stored
    logging.info(f"Logging to {log_file_path}")

def is_gzipped(file_path):
    """
    Check if a file is gzipped.
    """
    with open(file_path, 'rb') as f:
        return f.read(2) == b'\x1f\x8b'


def run_command(command):
    """Execute a shell command provided as a list, handling errors."""
    # logging.info("Executing command: %s", ' '.join(command))
    try:
        subprocess.run(command, check=True)
        logging.info("Command executed successfully: %s", ' '.join(command))
    except subprocess.CalledProcessError as e:
        logging.error("An error occurred while executing the command: %s", str(e))
        raise
        

def run_pints_visualizer(bam_file, output_prefix, mapq_threshold=255, experiment_type="CoPRO",
                         filters=["U13369", "chrM", "EBV", "_", "Mycoplasma"]):
    """
    Execute the pints_visualizer tool with specified parameters.

    Parameters:
    bam_file (str): Path to the input BAM file.
    output_prefix (str): Prefix for the output files.
    mapq_threshold (int, optional): Mapping quality threshold (default: 255).
    experiment_type (str, optional): Experiment type (default: "CoPRO").
    filters (list of str, optional): Filters for chromosome names (default: ["U13369", "chrM", "EBV", "_", "Mycoplasma"]).

    Example:
    run_pints_visualizer('sample.bam', 'output_prefix', mapq_threshold=30, experiment_type='ChIP-seq', filters=['chrM'])
    """

    # Construct the command to run pints_visualizer for TSS (5' signal)
    tss_command = [
        "pints_visualizer",
        "-b", bam_file,
        "-e", "CoPRO",  # Extracting the 5′ initiation signal (TSS)
        "--mapq-threshold", str(mapq_threshold),
        "--chromosome-start-with", "chr",
        "-o", output_prefix + "_5p",  # Output for 5' BigWig files
        "--filters"
    ] + filters

    # Construct the command to run pints_visualizer for pause sites (3' signal)
    pause_command = [
        "pints_visualizer",
        "-b", bam_file,
        "-e", "R1_5",  # Extracting the 3′ pause signal
        "-r",  # Reverse strand processing (ensuring correct read end extraction)
        "--mapq-threshold", str(mapq_threshold),
        "--chromosome-start-with", "chr",
        "-o", output_prefix + "_3p",  # Output for 3' BigWig files
        "--filters"
    ] + filters

    # Execute both commands
    run_command(tss_command)  # Run PINTS for 5′ end (TSS)
    run_command(pause_command)  # Run PINTS for 3′ end (pause site)

def quality_check(input_files, output_dir):
    """
    Perform a quality check with FastQC.
    Only runs FastQC if reports for all input files do not already exist.

    Parameters:
    - input_files (list of str): List of input FASTQ file paths.
    - output_dir (str): Directory where FastQC reports will be stored.

    Returns:
    - bool: True if FastQC was run or all reports already exist, False if an error occurred.
    """
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Determine which files need FastQC
    reports_to_process = []
    for input_file in input_files:
        base_name = os.path.basename(input_file).replace('.gz', '').replace('.fastq', '').rstrip('.fq')
        zip_report = os.path.join(output_dir, f"{base_name}_fastqc.zip")
        html_report = os.path.join(output_dir, f"{base_name}_fastqc.html")
        print(zip_report, html_report)
        if not (os.path.exists(zip_report) and os.path.exists(html_report)):
            reports_to_process.append(input_file)
    
    # Skip FastQC if all reports already exist
    if not reports_to_process:
        logging.info(f"Quality reports for all input files already exist in {output_dir}. Skipping FastQC.")
        return True

    # Run FastQC for missing reports
    logging.info(f"Running FastQC for files: {reports_to_process}")
    cmd = [
        'fastqc',  # FastQC command
        *reports_to_process,  # List of input files to process
        '-o', output_dir  # Output directory
    ]

    # Execute the FastQC command
    try:
        subprocess.run(cmd, check=True)
        logging.info("FastQC completed successfully.")
        return True
    except subprocess.CalledProcessError as e:
        logging.error(f"Error running FastQC: {e}")
        return False


def fastp_trim_umi_extraction(input_file1, input_file2, 
                              output_file1, output_file2,
                              json_report, html_report, 
                              umi_len, thread_count, adapter_seq1, adapter_seq2, 
                              overlap_len, length_required, trim_front1, trim_front2):
    """
    Reads trimming and UMI extraction using Fastp.
    Only performs trimming if the output files do not already exist.

    Parameters:
    input_file1 (str): Input FASTQ file for read 1.
    input_file2 (str): Input FASTQ file for read 2.
    output_file1 (str): Output FASTQ file for trimmed read 1.
    output_file2 (str): Output FASTQ file for trimmed read 2.
    json_report (str): Filename for the JSON report generated by Fastp.
    html_report (str): Filename for the HTML report generated by Fastp.
    umi_len (int): Length of UMI to extract.
    thread_count (int): Number of threads to use.
    adapter_seq1 (str): Adapter sequence for read 1.
    adapter_seq2 (str): Adapter sequence for read 2.
    overlap_len (int): Minimum overlap length required.
    length_required (int): Minimum length of reads required after trimming.
    trim_front1 (int): Number of bases to trim from the start of read 1.
    trim_front2 (int): Number of bases to trim from the start of read 2.
    """
    # Check if output files already exist
    if os.path.exists(output_file1) and os.path.exists(output_file2):
        logging.info(f"Trimmed files {output_file1} and {output_file2} already exist. Skipping trimming step.")
        return

    # Construct the command to run Fastp
    cmd = [
        'fastp',                           # The command to run the Fastp tool
        '-i', input_file1,                 # Input FASTQ file for read 1
        '-I', input_file2,                 # Input FASTQ file for read 2
        '-o', output_file1,                # Output FASTQ file for trimmed read 1
        '-O', output_file2,                # Output FASTQ file for trimmed read 2
        '-j', json_report,                 # Filename for the JSON report generated by Fastp
        '-h', html_report,                 # Filename for the HTML report generated by Fastp
        '-l', '18',                        # Only keep reads longer than 18 nt
        '--umi', '--umi_len={}'.format(umi_len), '--umi_loc=per_read',  # UMI extraction settings
        '--adapter_sequence', adapter_seq1,            # Adapter sequence for read 1
        '--adapter_sequence_r2', adapter_seq2,         # Adapter sequence for read 2
        '--overlap_len_require', str(overlap_len),     # Minimum overlap length required
        '--length_required', str(length_required),     # Minimum length of reads required after trimming
        '--trim_front1', str(trim_front1),             # Number of bases to trim from the start of read 1
        '--trim_front2', str(trim_front2),             # Number of bases to trim from the start of read 2
        '--low_complexity_filter', '-g', '-c', '-p',   # Filters and settings for low complexity reads
        '-w', str(thread_count)                        # Number of threads Fastp can use
    ]
    
    # Execute the constructed command
    run_command(cmd)

def download_latest_gencode_annotations(output_dir, file_format='gtf'):
    """
    Downloads and unzips the latest GENCODE annotations in GTF or GFF3 format with a progress bar.
    Only downloads if the file does not already exist in the output directory.

    Args:
    output_dir (str): Directory where the annotation files will be saved.
    file_format (str): Format of the annotation file to download ('gtf' or 'gff3'). Default is 'gtf'.

    Returns:
    str: Path to the uncompressed annotation file.

    Raises:
    ValueError: If an invalid file format is specified.
    ConnectionError: If the connection to the GENCODE FTP server fails.
    """
    # Validate file format input
    if file_format not in ['gtf', 'gff3']:
        raise ValueError("Invalid file format specified. Choose 'gtf' or 'gff3'.")

    # Base URL for GENCODE FTP server
    base_url = 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/'
    try:
        response = requests.get(base_url)
        response.raise_for_status()
    except requests.exceptions.HTTPError as err:
        raise ConnectionError(f"HTTP error occurred: {err}")
    except requests.exceptions.RequestException as err:
        raise ConnectionError(f"Error connecting to GENCODE FTP server: {err}")

    # Find the latest release number
    release_numbers = re.findall(r'release_([0-9]+)', response.text)
    if not release_numbers:
        raise ValueError("No release numbers found on the GENCODE page.")
    
    latest_release_number = max(map(int, release_numbers))

    # Construct file names and download URL
    file_suffix = 'gtf.gz' if file_format == 'gtf' else 'gff3.gz'
    compressed_file_name = f'gencode.v{latest_release_number}.primary_assembly.annotation.{file_suffix}'
    uncompressed_file_name = compressed_file_name[:-3]  # Remove .gz
    compressed_local_filename = os.path.join(output_dir, compressed_file_name)
    uncompressed_local_filename = os.path.join(output_dir, uncompressed_file_name)

    # Check if the uncompressed file already exists
    if os.path.exists(uncompressed_local_filename):
        logging.info(f"The file {uncompressed_local_filename} already exists. Skipping download.")
        return uncompressed_local_filename

    # Download the compressed file with a progress bar
    download_url = f'{base_url}release_{latest_release_number}/{compressed_file_name}'
    with requests.get(download_url, stream=True) as r:
        r.raise_for_status()
        total_size_in_bytes = int(r.headers.get('content-length', 0))
        chunk_size = 8192
        progress_bar = tqdm(total=total_size_in_bytes, unit='iB', unit_scale=True, desc="Downloading GENCODE annotations")
        with open(compressed_local_filename, 'wb') as f:
            for chunk in r.iter_content(chunk_size=chunk_size): 
                progress_bar.update(len(chunk))
                f.write(chunk)
        progress_bar.close()

    # Unzip the file and handle potential errors
    try:
        with gzip.open(compressed_local_filename, 'rb') as f_in:
            with open(uncompressed_local_filename, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
    except IOError as e:
        raise IOError(f"Failed to unzip file: {e}")

    # Optionally, remove the compressed file after extraction
    os.remove(compressed_local_filename)

    return uncompressed_local_filename

def download_genome(url, output_path):
    """
    Download a genome reference file from a specified URL with progress tracking,
    unzip it if necessary, and return the path to the uncompressed file.

    Parameters:
    url (str): URL to download the file from.
    output_path (str): Local path to save the downloaded file.

    Returns:
    str: The path to the uncompressed genome file.
    """
    # Configure logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    # Determine the output paths
    compressed_output_path = output_path + ".gz" if not output_path.endswith('.gz') else output_path
    uncompressed_output_path = output_path.rstrip('.gz')

    # Send a HTTP request to the server and get a response object with stream=True for large files
    response = requests.get(url, stream=True)
    response.raise_for_status()  # Raises an HTTPError for bad responses

    # Get the total file size from header
    total_size_in_bytes = int(response.headers.get('content-length', 0))

    # Open the output file and write chunks of data
    with open(compressed_output_path, 'wb') as f, tqdm(
        desc=f"Downloading {compressed_output_path}",
        total=total_size_in_bytes,
        unit='iB',
        unit_scale=True,
        unit_divisor=1024,
    ) as bar:
        for chunk in response.iter_content(chunk_size=8192):
            f.write(chunk)
            bar.update(len(chunk))

    # Log the completion of the download
    logging.info(f"Download completed: {compressed_output_path}")

    # Unzip the file
    with gzip.open(compressed_output_path, 'rb') as f_in, open(uncompressed_output_path, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

    # Remove the compressed file after extraction
    os.remove(compressed_output_path)

    # Log the completion of the unzipping
    logging.info(f"File uncompressed to: {uncompressed_output_path}")

    return uncompressed_output_path

def download_rDNA_curl(accession, output_path):
    """
    Download an rDNA sequence using curl from NCBI's website.

    Parameters:
    accession (str): The accession number of the rDNA sequence.
    output_path (str): Where to save the downloaded sequence.
    """
    url = f"https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id={accession}&db=nuccore&report=fasta"
    cmd = f"curl -o {output_path} \"{url}\""
    subprocess.run(cmd, shell=True)
    logging.info(f"Download completed: {output_path}")

def generate_star_index(genome_dir, fasta_files, gtf_file, run_threads):
    """
    Generate a STAR genome index.

    Parameters:
    star_path (str): Path to the STAR executable. (now removed)
    run_threads (int): Number of threads to use for indexing.
    genome_dir (str): Output directory for the genome indices.
    fasta_files (list of str): List of paths to genome FASTA files.
    gtf_file (str): Path to the GTF file.
    """
    logging.info(f"Starting STAR genome index generation in {genome_dir} with {run_threads} threads.")
    logging.info(f"Using FASTA files: {', '.join(fasta_files)}")
    logging.info(f"Using GTF file: {gtf_file}")
    
    # Construct the command to run STAR for genome index generation
    cmd = [
        'STAR',                             # The command to run the STAR tool
        '--runMode', 'genomeGenerate',      # Specify that STAR should run in genome generation mode
        '--runThreadN', str(run_threads),   # Number of threads to use for indexing
        '--genomeDir', genome_dir,          # Output directory for the genome indices
        '--genomeFastaFiles', *fasta_files, # List of paths to genome FASTA files
        '--sjdbGTFfile', gtf_file           # Path to the GTF file
    ]
    
    logging.info("Executing STAR command: " + ' '.join(cmd))
    
    # Execute the constructed command
    run_command(cmd)
    
    logging.info("STAR genome index generation completed successfully.")


def run_star_alignment(genome_dir, output_prefix,
                       read_files, 
                       align_mates_gap_max=1000,
                       align_intron_max=1000, # disabled by default
                       out_filter_multimap_nmax=10, 
                       out_filter_mismatch_nmax=1,
                       out_filter_multimap_score_range=0,
                       quant_mode='GeneCounts',
                       out_sam_attributes='All', 
                       genome_load='LoadAndKeep', 
                       limit_bam_sort_ram=30000000000,
                       run_threads = 8):
    """
    Align reads with STAR using specified parameters, only if the output files do not already exist.

    This function constructs and executes a command to run the STAR 
    (Spliced Transcripts Alignment to a Reference) tool for aligning 
    RNA-seq reads to a reference genome. The alignment parameters can 
    be customized using the function arguments.

    Parameters:
    run_threads (int): Number of threads to use for alignment.
    genome_dir (str): Directory containing the genome indices.
    output_prefix (str): Prefix for the output files.
    read_files (tuple of str): Input read files (paired-end reads), 
                               specified as a tuple of two file paths.
    align_mates_gap_max (int, optional): Maximum gap between paired-end 
                                         reads. Default is 1000.
    align_intron_max (int, optional): Maximum intron length. Default is 1000.
    out_filter_multimap_nmax (int, optional): Maximum number of multiple 
                                              alignments allowed. Default is 10.
    out_filter_mismatch_nmax (int, optional): Maximum number of mismatches 
                                              allowed. Default is 1.
    out_filter_multimap_score_range (int, optional): Score range for 
                                                     multi-mapping reads. 
                                                     Default is 0.
    quant_mode (str, optional): Quantification mode (e.g., 'GeneCounts'). 
                                Default is 'GeneCounts'.
    out_sam_attributes (str, optional): Additional SAM attributes to output. 
                                        Default is 'All'.
    genome_load (str, optional): Genome loading mode. Default is 'LoadAndKeep'.
    limit_bam_sort_ram (int, optional): Maximum RAM for BAM sorting. 
                                        Default is 30000000000 (30 GB).

    Notes:
    - Ensure that STAR is installed and accessible from the command line.
    - The input read files should be gzipped FASTQ files.
    """
    # Define the expected output BAM file
    expected_bam_file = f"{output_prefix}Aligned.sortedByCoord.out.bam"

    # Check if the output BAM file already exists
    if os.path.exists(expected_bam_file):
        logging.info(f"STAR alignment output {expected_bam_file} already exists. Skipping STAR alignment.")
        return

    # Construct the command to run STAR
    cmd = [
        'STAR',                                      # The command to run the STAR tool
        '--runThreadN', str(run_threads),            # Number of threads to use for alignment
        '--readFilesCommand', 'zcat',                # Command to decompress the input files (zcat for gzipped files)
        '--genomeDir', genome_dir,                   # Directory containing the genome indices
        '--readFilesIn', read_files[0], read_files[1], # Input read files (paired-end reads)
        '--outFileNamePrefix', output_prefix,        # Prefix for the output files
        '--alignMatesGapMax', str(align_mates_gap_max), # Maximum gap between paired-end reads
        '--alignIntronMax', str(align_intron_max),   # Maximum intron length
        '--outFilterMultimapNmax', str(out_filter_multimap_nmax), # Maximum number of multiple alignments allowed
        '--outFilterMismatchNmax', str(out_filter_mismatch_nmax), # Maximum number of mismatches allowed
        '--outFilterMultimapScoreRange', str(out_filter_multimap_score_range), # Score range for multi-mapping reads
        '--outSAMtype', 'BAM', 'SortedByCoordinate', # Output format and sorting order (BAM, sorted by coordinate)
        '--quantMode', quant_mode,                   # Quantification mode (e.g., GeneCounts)
        '--outSAMattributes', out_sam_attributes,    # Additional SAM attributes to output
        '--genomeLoad', genome_load,                 # Genome loading mode
        '--limitBAMsortRAM', str(limit_bam_sort_ram) # Maximum RAM for BAM sorting
    ]
    
    # Execute the constructed command
    logging.info(f"Running STAR alignment with command: {' '.join(cmd)}")
    run_command(cmd)

    # Log completion of STAR alignment
    if os.path.exists(expected_bam_file):
        logging.info(f"STAR alignment completed. Output saved to {expected_bam_file}.")
    else:
        logging.error(f"STAR alignment did not produce the expected output file: {expected_bam_file}.")

def index_bam(bam_file):
    """
    Index the BAM file using samtools.
    Only indexes the BAM file if the index file does not already exist.

    Parameters:
    bam_file (str): Path to the BAM file to be indexed.
    """
    # Construct the expected index file path (BAM index file typically has the same name with '.bai' appended)
    index_file = bam_file + '.bai'

    # Check if the index file already exists
    if os.path.exists(index_file):
        logging.info(f"BAM index file {index_file} already exists. Skipping indexing.")
        return

    # Construct the samtools index command
    cmd = ['samtools', 'index', bam_file]

    # Execute the indexing command
    try:
        logging.info(f"Indexing BAM file: {bam_file}")
        run_command(cmd)
    except subprocess.CalledProcessError as e:
        logging.error(f"An error occurred while indexing BAM file {bam_file}: {e}")
        raise

def filter_alignments(input_bam, output_bam, chromosomes=None, min_quality=255, threads=8):
    """
    Filters a BAM file by mapping quality and optionally by specific chromosomes.
    Only runs if the output BAM file does not already exist.

    Parameters:
    - input_bam (str): Path to the input BAM file.
    - output_bam (str): Path to the output filtered BAM file.
    - chromosomes (list of str, optional): List of chromosome names to keep (e.g., ['chr1', 'chr2', ..., 'chrX', 'chrY', 'chrM']).
                                           If None, all chromosomes are retained.
    - min_quality (int, optional): Minimum mapping quality (MAPQ) for filtering. Default: 255 (unique reads only).
    - threads (int, optional): Number of threads to use. Default: 8.

    Returns:
    - bool: True if filtering succeeds or file already exists, False if an error occurs.
    """
    if not shutil.which("samtools"):
        logging.error("Error: samtools is not installed or not in PATH. Install it with 'conda install -c bioconda samtools'")
        return False

    if not os.path.exists(input_bam):
        logging.error(f"Input BAM file not found: {input_bam}")
        return False

    # Skip processing if the output BAM file already exists
    if os.path.exists(output_bam) and os.path.exists(output_bam + ".bai"):
        logging.info(f"Filtered BAM file already exists: {output_bam}. Skipping filtering.")
        return True

    # Construct the samtools command for filtering by MAPQ
    cmd = [
        "samtools",
        "view", "-b",
        "-f", "3",  # Include only properly paired reads
        "-q", str(min_quality),  # Filter by mapping quality
        "-@", str(threads), input_bam  # Specify input BAM after options
    ]

    # If chromosome filtering is enabled, add chromosome inclusion filter
    if chromosomes:
        cmd += ["-o", output_bam] + chromosomes
    else:
        cmd += ["-o", output_bam]

    try:
        logging.info(f"Filtering {input_bam} (MAPQ ≥ {min_quality}, Chromosomes: {chromosomes if chromosomes else 'All'})...")
        subprocess.run(cmd, check=True)

        # Index the filtered BAM file
        index_cmd = ["samtools", "index", output_bam]
        subprocess.run(index_cmd, check=True)

        logging.info(f"Filtered BAM saved to: {output_bam}")
        return True
    except subprocess.CalledProcessError as e:
        logging.error(f"Error filtering BAM file: {e}")
        return False


def remove_pcr_duplicates(input_bam, output_bam, stats_prefix, umi_separator=':', paired=True):
    """
    Run UMI-tools deduplication on BAM files.
    Only performs deduplication if the output BAM file does not already exist.

    Parameters:
    input_bam (str): Path to the input BAM file.
    output_bam (str): Path to the output BAM file (deduplicated).
    stats_prefix (str): Prefix for the statistics file generated by UMI-tools.
    umi_separator (str): Separator used for UMIs in the BAM file.
    paired (bool): Indicates whether the reads are paired-end or single-end.
    """
    # Check if the deduplicated BAM file already exists
    if os.path.exists(output_bam):
        logging.info(f"Deduplicated BAM file {output_bam} already exists. Skipping UMI-tools deduplication.")
        return

    # Construct the UMI-tools deduplication command
    cmd = [
        'umi_tools', 'dedup',
        '--unpaired-reads=discard',
        '--umi-separator=' + umi_separator,
        '--paired' if paired else '--unpaired',
        '-I', input_bam,
        '--output-stats=' + stats_prefix,
        '-S', output_bam
    ]

    # Run the command
    try:
        logging.info("Running UMI-tools with command: %s", ' '.join(cmd))
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"An error occurred while running UMI-tools: {e}")
        raise

def generate_end_bams(input_bam, bam_read1, bam_read2, threads=4):
    """Generate BAM files for the read1 and read2 from paired-end reads."""

    '''
    Flag -f 66: This option filters reads based on the following criteria:
0x40 (64): First in a pair – the read is the first read in a paired-end sequence.
0x02 (2): Properly paired – as above, indicating proper alignment with its mate.
This selection typically captures the first read in a properly paired sequence, which, for certain library types like standard Illumina reads, will be the 3' end of the fragment if the library is prepared such that read 1 sequences the second strand.
    '''
    # Extract reads1 (left to right)
    cmd_read1 = [
        'samtools',          # samtools command
        'view',              # subcommand to view and filter BAM file
        '-f', '66',          # include only reads that are the first in a pair and properly aligned
        '--write-index',     # write an index file for the output BAM
        '-@', str(threads),  # number of threads to use
        '-o', bam_read1,     # output BAM file path for reads1
        input_bam            # input BAM file path
    ]

    logging.info("Extracting reads1...")
    run_command(cmd_read1)

    '''
    Flag -f 130: This option filters for reads that meet all of the following criteria (using binary bitwise AND with the flag 130):
0x80 (128): Second in a pair – the read is the second read in a paired-end sequence.
0x02 (2): Properly paired – the read is properly paired with its mate according to the aligner.
The combination results in reads that are the second in a pair and are properly paired. The reads must satisfy both conditions.
    '''
    # Extract reads2 (right to left)
    cmd_read2 = [
        'samtools',          # samtools command
        'view',              # subcommand to view and filter BAM file
        '-f', '130',         # include only reads that are the second in a pair and properly aligned
        '--write-index',     # write an index file for the output BAM
        '-@', str(threads),  # number of threads to use
        '-o', bam_read2,     # output BAM file path for reads2
        input_bam            # input BAM file path
    ]

    logging.info("Extracting reads2...")
    run_command(cmd_read2)

def fetch_chromosome_sizes(human_genome_path, genome_size_file):
    """
    Fetches the chromosome sizes from a human genome fasta file using pyfaidx and writes them to a file.

    Parameters:
    human_genome_path (str): Path to the human genome fasta file.
    genome_size_file (str): Output file path where chromosome sizes will be saved.
    """
    cmd = f'faidx {human_genome_path} -i chromsizes > {genome_size_file}'
    try:
        # Execute the command
        subprocess.run(cmd, shell=True, check=True)
        logging.info(f"Chromosome sizes were written to {genome_size_file}")
    except subprocess.CalledProcessError as e:
        logging.info(f"An error occurred: {e}")

def safe_file_move(source, destination):
    """
    Safely move a file from the source to the destination.

    This function attempts to copy a file from the source path to the destination path.
    If the copy is successful, it removes the original file from the source path.
    If any error occurs during the process, it logs an error message.

    Parameters:
    source (str): The path to the source file that needs to be moved.
    destination (str): The path to the destination where the file should be moved.

    Example:
    safe_file_move('/path/to/source/file.txt', '/path/to/destination/file.txt')

    The function will attempt to copy 'file.txt' from '/path/to/source/' to '/path/to/destination/'.
    If the copy is successful, 'file.txt' will be removed from '/path/to/source/'.
    If an error occurs, it will log the error message.

    Notes:
    - Ensure that the source file exists and the destination directory is writable.
    - This function uses `shutil.copy` for copying and `os.remove` for removing files.
    - It also logs the success or failure of the file move operation.
    """
    try:
        shutil.copy(source, destination)
        os.remove(source)
        logging.info(f"Successfully moved {source} to {destination}")
    except Exception as e:
        logging.error(f"Failed to move {source} to {destination}: {str(e)}")

def parse_bed_file(bed_file):
    """
    Parse a BED file to extract feature lengths for RPKM normalization using NumPy vectorization.
    
    Parameters:
    bed_file (str): Path to the BED file containing gene coordinates.
    
    Returns:
    dict: Dictionary mapping chromosome region keys (chrom, start, end) to feature lengths (bp).
    """
    data = np.loadtxt(bed_file, dtype=str, delimiter='\t', usecols=(0, 1, 2))
    chroms = data[:, 0]
    starts = data[:, 1].astype(int)
    ends = data[:, 2].astype(int)
    feature_lengths = {(chrom, start, end): end - start for chrom, start, end in zip(chroms, starts, ends)}
    return feature_lengths

def normalize_bedgraph(file, norm_type, is_minus_strand=False, features_file=None):
    """
    Normalize a BedGraph file by adjusting coverage values according to the specified normalization method.
    
    Parameters:
    file (str): Path to the input BedGraph file.
    norm_type (str): Normalization method. Options include:
        - 'CPM': Counts Per Million - normalizes by total read depth scaled to 1 million.
        - 'RPM': Reads Per Million - equivalent to CPM, used interchangeably.
        - 'BPM': Bins Per Million - normalizes by total read depth scaled to 1 billion.
        - 'RPKM': Reads Per Kilobase per Million - normalizes read depth considering feature length.
    is_minus_strand (bool): Whether the BedGraph represents the minus strand.
        - If True, ensures that normalized values remain negative.
        - If False, values remain positive.
    features_file(str, optional): Path to a BED file containing gene coordinates for RPKM normalization.
        - If None, RPKM normalization cannot be performed.
    
    Returns:
    str: Path to the normalized BedGraph file.
    
    The function reads the BedGraph file, normalizes the coverage column using the chosen method,
    and preserves negative values for the minus strand to maintain strand-specific signal polarity.
    """
    data = np.loadtxt(file, dtype=object, delimiter='\t')
    chroms = data[:, 0]
    starts = data[:, 1].astype(int)
    ends = data[:, 2].astype(int)
    coverage = data[:, 3].astype(float)
    
    # Compute total reads
    total_reads = np.sum(np.abs(coverage))  # Use absolute values to maintain negative values

    if total_reads == 0:
        logging.warning("Total reads are zero, skipping normalization.")
        return file

    # Load feature lengths from BED file if RPKM is requested
    feature_lengths = None
    if norm_type == 'RPKM':
        if features_file is None:
            raise ValueError("RPKM normalization requires a features file. Please provide a valid --features_file.")
        feature_lengths = parse_bed_file(features_file)

    # Apply Normalization
    if norm_type in ['CPM', 'RPM']:
        coverage = (coverage / total_reads) * 1e6 # Normalize to per million reads
    elif norm_type == 'BPM':
        coverage = (coverage / total_reads) * 1e9 # Normalize to per billion reads
    elif norm_type == 'RPKM' and feature_lengths is not None:
        # Vectorized lookup for region lengths
        region_keys = np.core.defchararray.add(np.core.defchararray.add(chroms, ':'), starts.astype(str))
        region_lengths = np.array([feature_lengths.get(key, end - start) for key, start, end in zip(region_keys, starts, ends)])
        
        # Avoid division by zero
        valid_regions = region_lengths > 0
        coverage[valid_regions] = (coverage[valid_regions] / (region_lengths[valid_regions] / 1000)) / (total_reads / 1e6)
    
    if is_minus_strand:
        coverage = -np.abs(coverage)  # Preserve negative values for minus strand
    
    data[:, 3] = coverage.astype(str)
    norm_file = f"{file}_{norm_type}.bedgraph"
    np.savetxt(norm_file, data, fmt='%s', delimiter='\t')
    return norm_file

def generate_bedgraph_and_bigwig(bam_file, output_prefix, genome_size_file, end='5', threads=4, normalize=False, normalization='CPM', features_file=None):
    """
    Generate BedGraph files from a BAM file and convert them to BigWig format for both strands.
    Only generates files if they do not already exist.

    Parameters:
    bam_file (str): Path to the input BAM file.
    output_prefix (str): Prefix for the output files.
    genome_size_file (str): Path to the genome size file required for BigWig conversion.
    threads (int): Number of threads to use for processing.
    end (str): '5' or '3', indicating the end for processing.
    normalize (bool): Whether to normalize the BedGraph file.
    normalization (str): Normalization method ('CPM', 'RPM', 'BPM', 'None').
    """
    # Define output file paths
    plus_bedgraph_file = f"{output_prefix}_plus.bedgraph"
    minus_bedgraph_file = f"{output_prefix}_minus.bedgraph"
    
    plus_bigwig_file = f"{output_prefix}_pl_{normalization}.bw" if normalize else f"{output_prefix}_pl.bw"
    minus_bigwig_file = f"{output_prefix}_mn_{normalization}.bw" if normalize else f"{output_prefix}_mn.bw"

    # Check if BigWig files already exist
    if os.path.exists(plus_bigwig_file) and os.path.exists(minus_bigwig_file):
        logging.info(f"BigWig files {plus_bigwig_file} and {minus_bigwig_file} already exist. Skipping generation.")
        return

    # Setup temporary file paths for bedGraph files
    plus_bedgraph = tempfile.NamedTemporaryFile(delete=False)
    minus_bedgraph = tempfile.NamedTemporaryFile(delete=False)

    # plus_tmp = tempfile.NamedTemporaryFile(delete=False)
    # minus_tmp = tempfile.NamedTemporaryFile(delete=False)
    
    # # Use the .name property (string paths)
    # plus_bedgraph = plus_tmp.name
    # minus_bedgraph = minus_tmp.name

    # Set BedTools commands based on read type
    if 'read2' in bam_file:  # Read 2 is initiation
        cmd_plus_bedgraph = ['bedtools', 'genomecov', '-ibam', bam_file, '-bg', '-strand', '+', '-5' if end == '5' else '-3']
        cmd_minus_bedgraph = ['bedtools', 'genomecov', '-ibam', bam_file, '-bg', '-strand', '-', '-5' if end == '5' else '-3']
    elif 'read1' in bam_file:
        cmd_plus_bedgraph = ['bedtools', 'genomecov', '-ibam', bam_file, '-bg', '-strand', '-', '-5' if end == '5' else '-3']
        cmd_minus_bedgraph = ['bedtools', 'genomecov', '-ibam', bam_file, '-bg', '-strand', '+', '-5' if end == '5' else '-3']

    # Generate BedGraph for plus strand
    logging.info("Generating BedGraph for plus strand...")
    with open(plus_bedgraph.name, 'w') as f:
        subprocess.run(cmd_plus_bedgraph, stdout=f, check=True)

    # Generate BedGraph for minus strand
    logging.info("Generating BedGraph for minus strand...")
    with open(minus_bedgraph.name, 'w') as f:
        subprocess.run(cmd_minus_bedgraph, stdout=f, check=True)

    # Convert minus strand values to negative using awk
    minus_bedgraph_temp = f"{minus_bedgraph.name}_neg"
    awk_cmd = f"awk '{{ $4 = -$4; print }}' {minus_bedgraph.name} > {minus_bedgraph_temp}"
    subprocess.run(awk_cmd, shell=True, check=True)
    os.replace(minus_bedgraph_temp, minus_bedgraph.name)

    # Normalize BedGraph files if required
    if normalize and normalization in ['CPM', 'RPM', 'BPM', 'RPKM']:
        logging.info(f"Normalizing BedGraph using {normalization}...")
        if normalization == 'RPKM' and not features_file:
            raise ValueError("RPKM normalization requires a features file. Please provide a valid --features_file.")

        plus_bedgraph = normalize_bedgraph(plus_bedgraph, normalization, is_minus_strand=False, features_file=features_file)
        minus_bedgraph = normalize_bedgraph(minus_bedgraph, normalization, is_minus_strand=True, features_file=features_file)
    else:
        # Convert tempfile objects to their paths if normalization is skipped
        plus_bedgraph = plus_bedgraph.name if hasattr(plus_bedgraph, "name") else plus_bedgraph
        minus_bedgraph = minus_bedgraph.name if hasattr(minus_bedgraph, "name") else minus_bedgraph

    # Move temporary BedGraph files to final locations
    safe_file_move(plus_bedgraph, plus_bedgraph_file)
    safe_file_move(minus_bedgraph, minus_bedgraph_file)

    # Convert BedGraph to BigWig
    cmd_bigwig_plus = ['bedGraphToBigWig', plus_bedgraph_file, genome_size_file, plus_bigwig_file]
    cmd_bigwig_minus = ['bedGraphToBigWig', minus_bedgraph_file, genome_size_file, minus_bigwig_file]

    logging.info(f"Converting plus strand BedGraph to BigWig: {plus_bigwig_file}")
    subprocess.run(cmd_bigwig_plus)

    logging.info(f"Converting minus strand BedGraph to BigWig: {minus_bigwig_file}")
    subprocess.run(cmd_bigwig_minus)

    # Optional cleanup of intermediate BedGraph files
    if os.path.exists(plus_bedgraph_file):
        os.remove(plus_bedgraph_file)
    if os.path.exists(minus_bedgraph_file):
        os.remove(minus_bedgraph_file)

def cleanup_files(file_list):
    """
    Remove specified files and folders from the filesystem.

    This function iterates over a list of file or folder paths and attempts to remove each.
    It logs a message indicating whether the file or folder was successfully removed or if it
    did not exist. If an error occurs during the removal process, it logs an error message.

    Parameters:
    file_list (list): A list of file or folder paths to be deleted.

    Example:
    cleanup_files(['/path/to/file1.txt', '/path/to/folder1'])

    The function will attempt to remove 'file1.txt' and 'folder1' from their respective paths.
    If a file or folder does not exist, it will log a message indicating no action was taken.
    If an error occurs during the removal process, it will log the error message.

    Notes:
    - Ensure that the files or folders exist and the program has permission to delete them.
    - This function uses `os.remove` to delete files, `shutil.rmtree` to remove directories,
      and `os.path.exists` to check for existence.
    - It logs the success or failure of each deletion operation.
    """
    for path in file_list:
        try:
            if os.path.exists(path):
                if os.path.isfile(path):
                    os.remove(path)
                    logging.info(f"Successfully removed file: {path}")
                elif os.path.isdir(path):
                    shutil.rmtree(path)
                    logging.info(f"Successfully removed folder: {path}")
            else:
                logging.info(f"No action taken, file or folder does not exist: {path}")
        except Exception as e:
            logging.error(f"Failed to remove {path}: {e}")

def count_reads_fastq(fastq_file):
    """
    Count the number of reads in a FASTQ file.

    This function counts the number of reads in a FASTQ file. If the file is gzipped 
    (i.e., ends with '.gz'), it uses `zcat` to decompress the file and then counts 
    the number of lines divided by 4 (since each read in a FASTQ file spans 4 lines).

    Parameters:
    fastq_file (str): The path to the FASTQ file.

    Returns:
    str: The number of reads in the FASTQ file.

    Example:
    read_count = count_reads_fastq('/path/to/file.fastq.gz')

    The function will return the read count for the specified FASTQ file and log the count.

    Notes:
    - Ensure that the FASTQ file exists and is readable.
    - This function assumes that the input file is properly formatted as a FASTQ file.
    - It uses shell commands (`zcat` or `cat` and `wc -l`) to count lines and calculate the read count.

    """
    cmd = f"zcat {fastq_file} | echo $((`wc -l`/4))" if fastq_file.endswith('.gz') else f"cat {fastq_file} | echo $((`wc -l`/4))"
    result = subprocess.run(cmd, shell=True, text=True, capture_output=True)
    read_count = result.stdout.strip()
    logging.info(f"Read count for {fastq_file}: {read_count}")
    return read_count

def count_reads_bam(bam_file):
    """
    Count the number of reads in a BAM file using samtools.

    This function counts the number of reads in a BAM file by using the `samtools view -c` 
    command, which returns the number of alignments in the BAM file.

    Parameters:
    bam_file (str): The path to the BAM file.

    Returns:
    str: The number of reads in the BAM file.

    Example:
    read_count = count_reads_bam('/path/to/file.bam')

    The function will return the read count for the specified BAM file and log the count.

    Notes:
    - Ensure that the BAM file exists and is readable.
    - This function requires that samtools is installed and accessible from the command line.
    - It uses the `subprocess` module to run the `samtools` command and capture its output.

    """
    cmd = ['samtools', 'view', '-c', bam_file]
    result = subprocess.run(cmd, capture_output=True, text=True)
    read_count = result.stdout.strip()
    logging.info(f"Read count for {bam_file}: {read_count}")
    return read_count

def compute_fragment_size(bam_file, output_histogram, threads=32):
    """
    Compute fragment size distribution for paired-end BAM files using deepTools.

    Parameters:
    - bam_file (str): Path to the input sorted BAM file.
    - output_histogram (str): Path to the output histogram image file (e.g., 'fragment_size.png').
    - threads (int, optional): Number of threads to use. Default is 32.

    Returns:
    - bool: True if the function completes successfully, False otherwise.
    """
    # Check if BAM file exists
    if not os.path.exists(bam_file):
        logging.error(f"BAM file not found: {bam_file}")
        return False

    # Check if output already exists
    if os.path.exists(output_histogram):
        logging.info(f"Fragment size histogram already exists: {output_histogram}. Skipping computation.")
        return True

    # Construct the command for bamPEFragmentSize
    cmd = [
        'bamPEFragmentSize',
        '--bamfiles', bam_file,
        '--histogram', output_histogram,
        '--numberOfProcessors', str(threads),
        '--plotFileFormat', 'png'  # Ensures output is in PNG format
    ]

    # Execute the command
    try:
        logging.info(f"Running bamPEFragmentSize on {bam_file}...")
        subprocess.run(cmd, check=True)
        logging.info(f"Fragment size histogram saved to {output_histogram}")
        return True
    except subprocess.CalledProcessError as e:
        logging.error(f"Error while running bamPEFragmentSize: {e}")
        return False
    

def sort_bam(input_bam, output_bam=None, threads=8):
    """
    Sorts a BAM file using samtools.

    Parameters:
    - input_bam (str): Path to the input BAM file.
    - output_bam (str, optional): Path to the sorted BAM file. Defaults to input_bam with "_sorted" suffix.
    - threads (int): Number of threads to use (default: 8).

    Returns:
    - str: Path to the sorted BAM file if successful, None otherwise.
    """
    if output_bam is None:
        output_bam = input_bam.replace(".bam", "_sorted.bam")

    cmd = ["samtools", "sort", "-@", str(threads), "-o", output_bam, input_bam]

    logging.info(f"Sorting BAM file: {input_bam} -> {output_bam}")
    
    try:
        subprocess.run(cmd, check=True)
        return output_bam
    except subprocess.CalledProcessError as e:
        logging.error(f"Error sorting BAM {input_bam}: {e}")
        return None


def parse_args():
    """Parse command line arguments."""
    """
    Parse command line arguments for processing paired-end FASTQ files to generate BigWig files.

    This function uses argparse to define and parse command line arguments required for the script.
    It includes options for input files, output directories, UMI extraction, adapter trimming,
    STAR alignment, and reference genome URLs.

    Returns:
    argparse.Namespace: An object containing the parsed command line arguments.

    Command Line Arguments:
    -f1, --fastq1 (str): Input FASTQ file for read 1 (required).
    -f2, --fastq2 (str): Input FASTQ file for read 2 (required).
    -o, --output_dir (str): Directory for output files (default: 'out_dir').
    -g, --genome_dir (str): Directory for genome files (default: 'genome_dir').
    -p, --prefix (str): Prefix for output file names (required).
    --spike_ins (flag): Flag to use multiple FASTA files for generating STAR index.
    -u, --umi_length (int): Length of UMI to extract (default: 6).
    -t, --threads (int): Number of threads to use (default: 48).
    -as1, --adapter_seq1 (str): Adapter sequence for read 1 (default: "TGGAATTCTCGGGTGCCAAGG").
    -as2, --adapter_seq2 (str): Adapter sequence for read 2 (default: "GATCGTCGGACTGTAGAACTCTGAAC").
    -ol, --overlap_length (int): Minimum overlap length required (default: 18).
    -lr, --length_required (int): Minimum length of reads required after trimming (default: 18).
    -tf1, --trim_front1 (int): Number of bases to trim from the start of read 1 (default: 0).
    -tf2, --trim_front2 (int): Number of bases to trim from the start of read 2 (default: 0).
    --human_genome_url (str): URL for the human genome (default: "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz").
    --drosophila_genome_url (str): URL for the Drosophila genome (default: "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/215/GCA_000001215.4_Release_6_plus_ISO1_MT/GCA_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz").
    --align_mates_gap_max (int): Maximum gap between two mates (default: 1000).
    --out_filter_multimap_nmax (int): Maximum number of multiple alignments allowed for a read (default: 10).
    --out_filter_mismatch_nmax (int): Maximum number of mismatches in the alignment (default: 1).
    --out_filter_multimap_score_range (int): Score range for filtering multiple alignments (default: 0).
    --align_intron_max (int): Maximum intron size allowed in the alignment (default: 1000).
    --quant_mode (str): Quantification mode (default: 'GeneCounts').
    --out_sam_attributes (str): Attributes to include in the output SAM (default: 'All').
    --genome_load (str): Mode of genome loading into memory (default: 'LoadAndKeep').
    --limit_bam_sort_ram (int): RAM limit for sorting BAM files (default: 30000000000).

    Example:
    args = parse_args()

    The function will parse the command line arguments and store them in the `args` object,
    which can then be accessed as needed within the script.

    Notes:
    - Ensure that required arguments are provided when running the script.
    - This function uses argparse to handle command line argument parsing.
    """
    parser = argparse.ArgumentParser(description="Process paired-end fastq files to generate BigWig files for various strands and priming.")
    
    # Required Input Files: Provide multiple files for replicates otherwise 1 file for each input
    parser.add_argument('-f1', '--fastq1', nargs='+', required=True, help='List of input FASTQ files for read 1 (must match -f2)')
    parser.add_argument('-f2', '--fastq2', nargs='+', required=True, help='List of input FASTQ files for read 2 (must match -f1)')
    parser.add_argument('-o', '--output_dir', default = 'out_dir', help='Prefix for output files')
    parser.add_argument('-g', '--genome_dir', default = 'genome_dir', help='Prefix for genome folder') #
    parser.add_argument('-p', '--prefix', required=True, help='Prefix for output file names')

    ## Determine which FASTA files to use for reference
    parser.add_argument('--spike_ins', action='store_true', help='Flag to use multiple FASTA files for generating STAR index')

    ## Trimming and UMI
    parser.add_argument('-u', '--umi_length', type=int, default=6, help='Length of UMI to extract')
    parser.add_argument('-t', '--threads', type=int, default=4, help='Number of threads to use')
    parser.add_argument('-as1', '--adapter_seq1', default="TGGAATTCTCGGGTGCCAAGG", help='Adapter sequence for read 1')
    parser.add_argument('-as2', '--adapter_seq2', default="GATCGTCGGACTGTAGAACTCTGAAC", help='Adapter sequence for read 2')
    parser.add_argument('-ol', '--overlap_length', type=int, default=18, help='Minimum overlap length required')
    parser.add_argument('-lr', '--length_required', type=int, default=18, help='Minimum length of reads required after trimming')
    parser.add_argument('-tf1', '--trim_front1', type=int, default=0, help='Number of bases to trim from the start of read 1')
    parser.add_argument('-tf2', '--trim_front2', type=int, default=0, help='Number of bases to trim from the start of read 2')

    ## STAR index
    parser.add_argument('--human_genome_url', default = "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz", help="URL for the human genome")
    parser.add_argument('--drosophila_genome_url', default = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/215/GCA_000001215.4_Release_6_plus_ISO1_MT/GCA_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz", help="URL for the Drosophila genome")

    ## STAR alignment
    parser.add_argument('-si', '--star_index_dir', default = 'star_index_dir', help='Location for STAR index folder') 
    parser.add_argument('--align_mates_gap_max', type=int, default=1000, help='Maximum gap between two mates.')
    parser.add_argument('--out_filter_multimap_nmax', type=int, default=10, help='Maximum number of multiple alignments allowed for a read.')
    parser.add_argument('--out_filter_mismatch_nmax', type=int, default=1, help='Maximum number of mismatches in the alignment.')
    parser.add_argument('--out_filter_multimap_score_range', type=int, default=0, help='Score range for filtering multiple alignments.')
    parser.add_argument('--align_intron_max', type=int, default=1000, help='Maximum intron size allowed in the alignment.')
    parser.add_argument('--quant_mode', type=str, default='GeneCounts', help='Quantification mode (typically GeneCounts).')
    parser.add_argument('--out_sam_attributes', type=str, default='All', help='Attributes to include in the output SAM.')
    parser.add_argument('--genome_load', type=str, default='LoadAndKeep', help='Mode of genome loading into memory.')
    parser.add_argument('--limit_bam_sort_ram', type=int, default=30000000000, help='RAM limit for sorting BAM files.')

    ## Generate bigwig using PINTS
    parser.add_argument('--run_pints_preprocess', action='store_true', help='Run pints visualizer after processing.')
    parser.add_argument('--normalize_bws', action='store_true', help='Save normalized reads, defaults to CPM.')
    
    ## Bam to bigwig conversion
    parser.add_argument('-n', '--norm_type', choices=['CPM', 'RPM', 'BPM', 'None'], default='CPM',
                        help="Normalization method. Options: 'CPM', 'RPKM', 'BPM', or 'None'. Default: 'CPM'.")
    parser.add_argument('--bin_size', type=int, default=1, help='Bin size for smoothing. Default: 1 (base-pair resolution).')
    parser.add_argument("--features_file", type=str, default=None, help="Path to features file for RPKM normalization (required for RPKM method).")
    
    return parser.parse_args()

def main():
    args = parse_args()

    # Setup output directory
    output_dir = os.path.join(args.output_dir, args.prefix)
    setup_logging(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    ### STEP 1: PREPARE GENOME FILES ###
    if args.genome_dir != 'genome_dir':
        genome_dir = args.genome_dir
    else:
        genome_dir = os.path.join(args.output_dir, args.genome_dir) 
    os.makedirs(genome_dir, exist_ok=True)

    human_genome_path = os.path.join(genome_dir, "GRCh38_reference.fna")
    drosophila_genome_path = os.path.join(genome_dir, "dm6_genomic.fna")
    rdna_path = os.path.join(genome_dir, "rDNA_U13369.1.fa")

    # Check and download genome files
    if not os.path.exists(human_genome_path):
        download_genome(args.human_genome_url, human_genome_path)

    if not os.path.exists(rdna_path):
        download_rDNA_curl("U13369.1", rdna_path)

    # Generate chromosome sizes
    genome_size_file = os.path.join(genome_dir, 'hg38.chrom.sizes')
    if not os.path.exists(genome_size_file):
        fetch_chromosome_sizes(human_genome_path, genome_size_file)

    # Download the latest Gencode annotations
    gtf_path = download_latest_gencode_annotations(genome_dir, file_format='gtf')

    # Handle optional spike-ins
    if args.spike_ins:
        if not os.path.exists(drosophila_genome_path):
            download_genome(args.drosophila_genome_url, drosophila_genome_path)
        fastas = [human_genome_path, rdna_path, drosophila_genome_path]
    else:
        fastas = [human_genome_path, rdna_path]

    # Generate STAR index if not already present
    if args.star_index_dir != "star_index_dir":
        star_index_dir = args.star_index_dir
    else:
        star_index_dir = os.path.join(args.output_dir, 'star_index_dir') 
        if not os.path.exists(star_index_dir):
            os.makedirs(star_index_dir, exist_ok=True) 
            generate_star_index(star_index_dir, fastas, gtf_path, args.threads)

    ### STEP 2: READ INPUT FASTQ FILES ###
    if len(args.fastq1) != len(args.fastq2):
        raise ValueError("Error: The number of files in -f1 and -f2 must match.")
    fastq_pairs = list(zip(args.fastq1, args.fastq2))
    logging.info(f"Processing {len(fastq_pairs)} replicate pairs.")

    # Determine if multiple replicates exist
    single_replicate = len(fastq_pairs) == 1
    processed_bams = []  # Store deduplicated BAMs for merging

    ### STEP 3: PROCESS EACH REPLICATE ###
    for i, (fastq1, fastq2) in enumerate(fastq_pairs, 1):
        # rsample_prefix = f"{args.prefix}_rep{i}"
        sample_prefix = args.prefix if single_replicate else f"{args.prefix}_rep{i}"
        sample_output_dir = os.path.join(output_dir, sample_prefix)
        os.makedirs(sample_output_dir, exist_ok=True)

        logging.info(f"Processing replicate {i}: {fastq1}, {fastq2}")

        # Validate FASTQ files
        for file in [fastq1, fastq2]:
            if not is_gzipped(file):
                raise ValueError(f"Error: {file} is not a gzipped file.")

        # Perform quality check
        quality_check([fastq1, fastq2], sample_output_dir)
        count_reads_fastq(fastq1)
        count_reads_fastq(fastq2)

        # Trim and extract UMIs and generate quality reports
        trimmed_fastq1 = os.path.join(sample_output_dir, f"{sample_prefix}_trimmed_1.fq.gz")
        trimmed_fastq2 = os.path.join(sample_output_dir, f"{sample_prefix}_trimmed_2.fq.gz")
        fastp_trim_umi_extraction(
            fastq1, fastq2, trimmed_fastq1, trimmed_fastq2,
            os.path.join(sample_output_dir, f"{sample_prefix}.json"),
            os.path.join(sample_output_dir, f"{sample_prefix}.html"),
            args.umi_length,
            args.threads, 
            args.adapter_seq1, args.adapter_seq2,
            args.overlap_length, 
            args.length_required, 
            args.trim_front1, args.trim_front2
        )
        # Count reads in trimmed files
        count_reads_fastq(trimmed_fastq1)
        count_reads_fastq(trimmed_fastq2)
        
#         ## Perform quality check on trimmed files
#         quality_check([trimmed_fastq1, trimmed_fastq2], sample_output_dir)
        
        # Align reads with STAR
        star_out_dir = os.path.join(output_dir, args.prefix+'_star_output')
        if not os.path.exists(star_out_dir):
            os.makedirs(star_out_dir, exist_ok=True) 
        aligned_bam = os.path.join(star_out_dir, f"{sample_prefix}_Aligned.sortedByCoord.out.bam")

        run_star_alignment(
                    genome_dir=star_index_dir,
                    output_prefix=os.path.join(star_out_dir, sample_prefix + '_'),
                    read_files= (trimmed_fastq1, trimmed_fastq2), 
                    align_mates_gap_max=args.align_mates_gap_max,
                    out_filter_multimap_nmax=args.out_filter_multimap_nmax,
                    out_filter_mismatch_nmax=args.out_filter_mismatch_nmax,
                    out_filter_multimap_score_range=args.out_filter_multimap_score_range,
                    align_intron_max=args.align_intron_max,
                    quant_mode=args.quant_mode,
                    out_sam_attributes=args.out_sam_attributes,
                    genome_load=args.genome_load,
                    limit_bam_sort_ram=args.limit_bam_sort_ram,
                    run_threads=args.threads
                )
        index_bam(aligned_bam)
        count_reads_bam(aligned_bam)

        # Filter unique reads by mapping quality
        quality_filtered_bam = os.path.join(star_out_dir, f"{sample_prefix}_qfiltered.bam")
        filter_alignments(aligned_bam, quality_filtered_bam, min_quality=255, threads=args.threads)
        count_reads_bam(quality_filtered_bam)
        index_bam(quality_filtered_bam)
        
        # Compute Fragment Size Distribution before filtering
        fragment_size_plot = os.path.join(sample_output_dir, f"{sample_prefix}_fragment_size.png")
        compute_fragment_size(aligned_bam, fragment_size_plot)

        ## Remove PCR duplicates
        dedup_bam = os.path.join(star_out_dir, sample_prefix+'_dedup.bam')    
        umi_separator = ":" 
        stats_file_prefix = os.path.join(sample_output_dir, sample_prefix)
        remove_pcr_duplicates(quality_filtered_bam, dedup_bam, stats_file_prefix, umi_separator)
        count_reads_bam(dedup_bam)
        index_bam(dedup_bam)
        
        # Filter for chromosomes of interest
        dedup_filtered_bam = os.path.join(star_out_dir, sample_prefix+'_dedup_chromfiltered.bam')
        filter_alignments(dedup_bam, dedup_filtered_bam, chromosomes=["chr" + str(i) for i in range(1, 23)] + ["chrX", "chrY", "chrM"])
        count_reads_bam(dedup_filtered_bam)
        index_bam(dedup_filtered_bam)

        # Sort bam file
        dedup_sorted_bam = os.path.join(star_out_dir, sample_prefix+'_dedup_chromfiltered.bam')
        sort_bam(dedup_filtered_bam, dedup_sorted_bam, threads = args.threads)
        
        # Append each dedup bam to list 
        processed_bams.append(dedup_sorted_bam)

    ### STEP 4: MERGE REPLICATES IF MORE THAN 1 REPLICATES ###
    if len(processed_bams) > 1:
        merged_bam = os.path.join(star_out_dir, f"{args.prefix}_merged.bam")
        merge_command = ["samtools", "merge", "-f", "-@", str(args.threads), merged_bam] + processed_bams
        logging.info(f"Merging {len(processed_bams)} replicates into {merged_bam}...")
        subprocess.run(merge_command, check=True)
        index_bam(merged_bam)
        final_bam = merged_bam
    else:
        logging.info("Only one replicate provided, skipping merging.")
        final_bam = processed_bams[0] if processed_bams else None  # Use single replicate BAM if available

    # Run PINTS visualizer to generate bigwig files
    if args.run_pints_preprocess:
        run_pints_visualizer(final_bam, os.path.join(output_dir,args.prefix+"_pints_output"))
        
    ### STEP 5: GENERATE BIGWIG FILES ###
    logging.info(f"Generating BigWig files from deduplicated BAM: {final_bam}")
    bam_read1 = final_bam.replace(".bam", "_read1.bam")
    bam_read2 = final_bam.replace(".bam", "_read2.bam")
    generate_end_bams(final_bam, bam_read1, bam_read2, threads=args.threads)

    ## Generate BigWig files (second read is initiation)
    features_file=args.features_file if args.norm_type == 'RPKM' else None
    generate_bedgraph_and_bigwig(bam_read2, os.path.join(output_dir, args.prefix), genome_size_file, end='5', threads=args.threads, normalize=args.normalize_bws, normalization=args.norm_type, features_file=features_file)
    generate_bedgraph_and_bigwig(bam_read1, os.path.join(output_dir, args.prefix + '_3p'), genome_size_file, end='5', threads=args.threads, normalize=args.normalize_bws, normalization=args.norm_type, features_file=features_file)

    # # Call cleanup with intermediate files
    # prefix = args.prefix
    # # Define intermediate files to remove
    # intermediate_files = [
    #     '_STARtmp',
    #     os.path.join(output_dir, sample_prefix, f"{sample_prefix}_trimmed_1.fq.gz"),
    #     # os.path.join(sample_output_dir, f"{sample_prefix}.json")
    #     os.path.join(output_dir, sample_prefix, f"{sample_prefix}_trimmed_2.fq.gz"),
    #     os.path.join(output_dir, sample_prefix, f"{sample_prefix}.json"),
    #     os.path.join(output_dir, sample_prefix, f"{sample_prefix}.html"),
    #     os.path.join(star_out_dir, f"{sample_prefix}_Aligned.sortedByCoord.out.bam"),
    #     os.path.join(star_out_dir, f"{sample_prefix}_Aligned.sortedByCoord.out.bam.bai"),
    #     os.path.join(star_out_dir, f"{sample_prefix}_qfiltered.bam"),
    #     os.path.join(star_out_dir, f"{sample_prefix}_qfiltered.bam.bai"),
    #     os.path.join(star_out_dir, f"{sample_prefix}_dedup.bam"),
    #     os.path.join(star_out_dir, f"{sample_prefix}_dedup.bam.bai"),
    #     os.path.join(star_out_dir, f"{sample_prefix}_dedup_chromfiltered.bam"),
    #     os.path.join(star_out_dir, f"{sample_prefix}_dedup_chromfiltered.bam.bai"),
    #     os.path.join(output_dir, f"{prefix}_plus.bedgraph"),
    #     os.path.join(output_dir, f"{prefix}_minus.bedgraph"),
    #     os.path.join(output_dir, f"{prefix}_3p_plus.bedgraph"),
    #     os.path.join(output_dir, f"{prefix}_3p_minus.bedgraph")
    # ]
    # # Perform cleanup
    # # cleanup_files(intermediate_files)
    
    logging.info("Processing complete!")

if __name__ == "__main__":
    main()


