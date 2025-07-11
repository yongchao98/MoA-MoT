#!/usr/bin/env python

import argparse
import pandas as pd

def filter_contaminants(input_file, output_file, log2fc_thresh, padj_thresh):
    """
    Filters a differential expression file to remove likely contaminants.

    Contaminants are defined as genes with a log2FoldChange < log2fc_thresh
    and padj < padj_thresh. This is based on the assumption that contamination
    is present in the reference/control group of the comparison.

    Args:
        input_file (str): Path to the input differential expression CSV/TSV file.
        output_file (str): Path to save the filtered CSV file.
        log2fc_thresh (float): The negative log2 fold change threshold.
        padj_thresh (float): The adjusted p-value threshold for significance.
    """
    # Detect separator (CSV or TSV)
    sep = ',' if input_file.endswith('.csv') else '\t'
    
    # Read the differential expression data
    try:
        de_data = pd.read_csv(input_file, sep=sep)
    except FileNotFoundError:
        print(f"Error: Input file not found at '{input_file}'")
        return
    except Exception as e:
        print(f"Error reading file: {e}")
        return

    # Ensure required columns are present
    required_cols = ['log2FoldChange', 'padj']
    if not all(col in de_data.columns for col in required_cols):
        print(f"Error: Input file must contain the columns: {', '.join(required_cols)}")
        return

    initial_gene_count = len(de_data)

    # Define the filtering condition for contaminants
    # Contaminants have high expression in the control group, so their
    # log2FC is significantly negative.
    is_contaminant = (de_data['log2FoldChange'] < log2fc_thresh) & (de_data['padj'] < padj_thresh)
    
    # Separate the contaminants and the clean data
    contaminant_genes = de_data[is_contaminant]
    filtered_data = de_data[~is_contaminant]

    contaminant_count = len(contaminant_genes)
    final_gene_count = len(filtered_data)
    
    # Print the summary of the filtering operation
    print("--- Contamination Filtering Summary ---")
    print("\nFiltering Rule (gene is a contaminant if):")
    print(f"log2FoldChange < {log2fc_thresh}")
    print("AND")
    print(f"padj < {padj_thresh}\n")
    
    print(f"Total genes in input file: {initial_gene_count}")
    print(f"Number of contaminating genes removed: {contaminant_count}")
    print(f"Number of genes remaining after filtering: {final_gene_count}\n")

    # Save the filtered data to the output file
    filtered_data.to_csv(output_file, index=False)
    print(f"Filtered data has been saved to: {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Filter contaminating genes from a differential expression results file. "
                    "This script assumes contaminants are overly represented in the control condition, "
                    "resulting in a large negative log2FoldChange."
    )
    parser.add_argument(
        "-i", "--input", 
        required=True, 
        help="Path to the input differential expression file (CSV or TSV). Must contain 'log2FoldChange' and 'padj' columns."
    )
    parser.add_argument(
        "-o", "--output", 
        required=True, 
        help="Path for the output filtered data file (CSV format)."
    )
    parser.add_argument(
        "--log2fc_threshold", 
        type=float, 
        default=-2.0,
        help="Log2 fold change threshold to identify contaminants. Genes below this value are considered. Default: -2.0"
    )
    parser.add_argument(
        "--padj_threshold", 
        type=float, 
        default=0.05,
        help="Adjusted p-value threshold for significance. Default: 0.05"
    )

    args = parser.parse_args()

    filter_contaminants(args.input, args.output, args.log2fc_threshold, args.padj_threshold)
    
    # Example of how to run from the command line:
    # python your_script_name.py --input /path/to/de_results.csv --output /path/to/filtered_results.csv