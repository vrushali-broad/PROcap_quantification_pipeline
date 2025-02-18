"""
Genome-wide Correlation Analysis of BigWig Files

Author: Vrushali D. Fangal  
Date: 2024-05-23  
Organization: YuLab, Cornell University  

Description:
This script calculates and visualizes genome-wide correlation between two BigWig files.
It bins the genomic data, extracts values, and computes a correlation scatter plot.
Multi-processing is used to speed up the per-chromosome binning process.

Features:
- Reads chromosome sizes from an input file.
- Extracts signal values from BigWig files in non-overlapping bins.
- Computes Pearson correlation between log-transformed values.
- Generates a scatter plot with a 45-degree reference line.
- Utilizes multiprocessing to speed up processing.

Usage:
```bash
python procap_correlation.py -f1 file1.bw -f2 file2.bw --bin_size 1000 --chrom_sizes_file chrom.sizes
```
"""

import pyBigWig
import numpy as np
import matplotlib.pyplot as plt
import argparse
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
import os
from scipy.stats import pearsonr
import logging

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

def load_chromosome_sizes(file_path):
    """Load chromosome sizes from a file."""
    logging.info("Loading chromosome sizes from file.")
    standard_chromosomes = {f"chr{i}" for i in list(range(1, 23)) + ["X", "Y"]}
    chrom_sizes = {}

    try:
        with open(file_path, "r") as file:
            for line in file:
                chrom_name, size = line.strip().split()[:2]
                if chrom_name in standard_chromosomes:
                    chrom_sizes[chrom_name] = int(size)
        logging.info("Chromosome sizes loaded successfully.")
    except Exception as e:
        logging.error(f"Error loading chromosome sizes: {e}")
        raise

    return chrom_sizes

def process_chromosome(args):
    # print(f"Process {os.getpid()} is handling {chrom}")  # Debugging print
    bw1_path, bw2_path, chrom, length, bin_size = args
    bins = []
    
    try:
        with pyBigWig.open(bw1_path) as bw1, pyBigWig.open(bw2_path) as bw2:
            if chrom not in bw1.chroms() or chrom not in bw2.chroms():
                return chrom, []  # Skip missing chromosomes

            num_bins = length // bin_size

            for i in range(num_bins):
                start = i * bin_size
                end = min(start + bin_size, length)
                bin_val1 = np.abs(np.nansum(bw1.values(chrom, start, end)))
                bin_val2 = np.abs(np.nansum(bw2.values(chrom, start, end)))
                if bin_val1 == 0 or bin_val2 == 0:
                    continue # Skip this window

                bins.append((bin_val1, bin_val2))

    except Exception as e:
        logging.error(f"Error processing chromosome {chrom}: {e}")
        return chrom, []

    return chrom, bins

def values_per_bin(bw1_path, bw2_path, chroms, bin_size):
    """Extract values per bin from both BigWig files for the exact same genomic regions."""
    bw1_name = os.path.basename(bw1_path)
    bw2_name = os.path.basename(bw2_path)

    all_data = []
    task_args = [(bw1_path, bw2_path, chrom, length, bin_size) for chrom, length in chroms.items()]

    try:
        with Pool(processes=cpu_count()) as pool:
            results = list(tqdm(pool.imap(process_chromosome, task_args), total=len(task_args), desc=f"Processing {bw1_name} vs {bw2_name}"))

            for _, bin_values in results:
                all_data.extend(bin_values)

    except Exception as e:
        logging.error(f"Error during multiprocessing: {e}")
        raise

    return zip(*all_data)

def plot_genome_correlation(bw1_path, bw2_path, chroms, bin_size):
    """Compute correlation and plot scatter plot for two BigWig files."""
    logging.info("Starting correlation plot generation.")

    try:
        data1, data2 = values_per_bin(bw1_path, bw2_path, chroms, bin_size)
    except Exception as e:
        logging.error(f"Error extracting values per bin: {e}")
        raise

    try:
        data1 = np.array(data1)
        data2 = np.array(data2)

        # Remove NaN and zero values
        valid_indices = (~np.isnan(data1)) & (~np.isnan(data2)) & (data1 > 0) & (data2 > 0)
        filtered_data1 = data1[valid_indices]
        filtered_data2 = data2[valid_indices]

        if len(filtered_data1) == 0 or len(filtered_data2) == 0:
            logging.error("No valid data points found for correlation.")
            return

        # Log-transform values
        log_data1 = np.log10(filtered_data1)
        log_data2 = np.log10(filtered_data2)

        # Compute Pearson correlation
        corr, _ = pearsonr(log_data1, log_data2)

        # Extract file names
        bw1_name = os.path.basename(bw1_path)
        bw2_name = os.path.basename(bw2_path)

        # Plot correlation scatter plot
        plt.figure(figsize=(5, 5))
        plt.scatter(log_data1, log_data2, alpha=0.1, s=2)
        plt.title(f"Genome-wide Correlation: {bw1_name} vs {bw2_name}", fontsize=10)
        plt.xlabel(f"Log10({bw1_name} counts)", fontsize=8)
        plt.ylabel(f"Log10({bw2_name} counts)", fontsize=8)

        # Set plot limits
        min_val, max_val = min(log_data1.min(), log_data2.min()), max(log_data1.max(), log_data2.max())
        padding = (max_val - min_val) * 0.05
        axis_limits = [min_val - padding, max_val + padding]

        plt.xlim(axis_limits)
        plt.ylim(axis_limits)
        plt.axis("equal")  # Ensures a 1:1 aspect ratio
        plt.plot(axis_limits, axis_limits, "r--", linewidth=1)  # 45-degree reference line
        plt.grid(True)

        plt.figtext(0.5, 0.99, f"Pearson Correlation: {corr:.2f}", ha="center", fontsize=12,
                    bbox={"facecolor": "grey", "alpha": 0.1, "pad": 5})

        plt.savefig(f"correlation_binsize_{bin_size}.png", bbox_inches="tight")
        logging.info("Correlation plot saved successfully.")

    except Exception as e:
        logging.error(f"Error during plotting: {e}")
        raise

def main():
    parser = argparse.ArgumentParser(description="Plot genome-wide correlation between two BigWig files.")
    parser.add_argument("-f1", "--bw1_path", type=str, required=True, help="Path to the first BigWig file.")
    parser.add_argument("-f2", "--bw2_path", type=str, required=True, help="Path to the second BigWig file.")
    parser.add_argument("--bin_size", type=int, default=1000, help="Size of bins for averaging.")
    parser.add_argument("--chrom_sizes_file", type=str, required=True, help="File containing chromosome sizes.")

    args = parser.parse_args()

    try:
        chrom_sizes = load_chromosome_sizes(args.chrom_sizes_file)
        plot_genome_correlation(args.bw1_path, args.bw2_path, chrom_sizes, args.bin_size)
    except Exception as e:
        logging.error(f"Error in main: {e}")
        raise

if __name__ == "__main__":
    main()

    
    
