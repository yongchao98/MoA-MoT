#!/usr/bin/env python
import pandas as pd
import sys
import argparse

def filter_contaminants(input_stream, output_stream, log2fc_threshold):
    """
    Filters a differential expression file to remove potential contaminants.

    Contaminants in the control group appear as genes with a large negative
    log2FoldChange. This script removes all genes below the specified
    threshold.

    Args:
        input_stream: A file-like object to read data from (e.g., sys.stdin).
        output_stream: A file-like object to write data to (e.g., sys.stdout).
        log2fc_threshold (float): The minimum log2FoldChange to keep.
    """
    try:
        # Detect separator by reading the first line
        first_line = input_stream.readline()
        if '\t' in first_line:
            sep = '\t'
        else:
            sep = ','
        # Reset stream and read the full file
        input_stream.seek(0)
        df = pd.read_csv(input_stream, sep=sep)

        # Ensure the log2FoldChange column exists
        if 'log2FoldChange' not in df.columns:
            print("Error: 'log2FoldChange' column not found in the input.", file=sys.stderr)
            sys.exit(1)

        initial_count = len(df)

        # The filtering condition is the core of the strategy
        # final_equation: log2FoldChange >= log2fc_threshold
        filtered_df = df[df['log2FoldChange'] >= log2fc_threshold].copy()

        final_count = len(filtered_df)
        removed_count = initial_count - final_count

        # Print an informational message to stderr (console)
        # This keeps the data output (stdout) clean for piping
        print(f"Filtering based on the equation: log2FoldChange >= {log2fc_threshold}", file=sys.stderr)
        print(f"Initial number of genes: {initial_count}", file=sys.stderr)
        print(f"Number of genes removed: {removed_count}", file=sys.stderr)
        print(f"Final number of genes: {final_count}", file=sys.stderr)

        # Print the filtered dataframe to standard output as a CSV
        filtered_df.to_csv(output_stream, index=False)

    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Filter contaminant genes from a differential expression file based on log2FoldChange.",
        epilog="Example Usage: cat dge_results.csv | python this_script.py > filtered_dge_results.csv"
    )
    parser.add_argument(
        '--threshold',
        type=float,
        default=-5.0,
        help='The minimum log2FoldChange value to keep. Default is -5.0.'
    )
    args = parser.parse_args()

    # The script reads from standard input and writes to standard output
    filter_contaminants(sys.stdin, sys.stdout, args.threshold)