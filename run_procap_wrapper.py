import os
import re
import subprocess
import argparse
import logging

# Set up logging
logging.basicConfig(
    filename="procap_wrapper.log",  # Log file name
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)

def find_fastq_pairs(fastq_dir):
    """Find pairs of FASTQ files in the given directory."""
    fastq_pattern = re.compile(r"(.+)_([12])\.fq\.gz$")
    fastq_files = sorted(f for f in os.listdir(fastq_dir) if f.endswith(".fq.gz"))

    file_pairs = {}
    for file in fastq_files:
        match = fastq_pattern.match(file)
        if match:
            sample_name, read_number = match.groups()
            if sample_name not in file_pairs:
                file_pairs[sample_name] = {}
            file_pairs[sample_name][read_number] = os.path.join(fastq_dir, file)

    return file_pairs

def run_procap_pipeline(fastq1, fastq2, output_dir, genome_dir, sample, args):
    """Run the PROcap Quantification Pipeline for a pair of FASTQ files."""
    # sample_output = os.path.join(output_dir, sample)
    # os.makedirs(sample_output, exist_ok=True)

    cmd = [
        "python", "PROcap_processing_pipeline.py",
        "-f1", fastq1, "-f2", fastq2,
        "-o", output_dir,
        "-p", sample,
        "-g", genome_dir,
        "--normalize_bws",
        "-t", str(args.threads)  # Pass thread count
    ]

    logging.info(f"Processing sample: {sample}")
    logging.info(f"FASTQ1: {fastq1}")
    logging.info(f"FASTQ2: {fastq2}")
    logging.info(f"Using {args.threads} threads")

    print(f"Running PROcap pipeline for: {sample}")
    subprocess.run(cmd, check=True)

def main():
    """Main function to parse arguments and process FASTQ pairs."""
    parser = argparse.ArgumentParser(description="Wrapper script for running PROcap Quantification Pipeline on paired FASTQ files.")
    parser.add_argument("-d", "--fastq_dir", required=True, help="Directory containing FASTQ files")
    parser.add_argument("-o", "--output_dir", default = 'Results',  help="Directory to store output files")
    parser.add_argument("-g", "--genome_dir", default = "genome_dir", help="Directory name for genome reference files")
    parser.add_argument("-t", "--threads", type=int, default=48, help="Number of threads to use (default: 48)")   

    args = parser.parse_args()

    # Ensure output directory exists
    os.makedirs(args.output_dir, exist_ok=True)

    # Find FASTQ pairs
    fastq_pairs = find_fastq_pairs(args.fastq_dir)

    # Process each FASTQ pair
    for sample, reads in fastq_pairs.items():
        if "1" in reads and "2" in reads:
            run_procap_pipeline(reads["1"], reads["2"], args.output_dir, args.genome_dir, sample, args)
        else:
            logging.warning(f"Skipping {sample}: Missing paired FASTQ files")
            print(f"Skipping {sample}: Missing paired FASTQ files")

if __name__ == "__main__":
    main()


