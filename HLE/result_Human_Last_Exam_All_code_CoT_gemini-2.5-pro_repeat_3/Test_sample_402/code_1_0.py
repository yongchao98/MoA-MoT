import pandas as pd
import argparse
import sys
import os

def filter_contaminants_by_logfc(file_path, log2fc_col, threshold):
    """
    Filters a differential expression results file to remove potential contaminants
    based on a log2 fold change threshold.

    Args:
        file_path (str): Path to the input CSV or TSV file.
        log2fc_col (str): The name of the column containing log2 fold change values.
        threshold (float): The log2FC value above which genes will be removed.
    """
    # --- 1. Input Validation and Data Loading ---
    if not os.path.exists(file_path):
        print(f"Error: Input file not found at '{file_path}'", file=sys.stderr)
        sys.exit(1)

    # Determine file separator based on extension
    if file_path.endswith('.csv'):
        separator = ','
    elif file_path.endswith('.tsv') or file_path.endswith('.txt'):
        separator = '\t'
    else:
        print(f"Error: Unsupported file format for '{file_path}'. Please use .csv or .tsv.", file=sys.stderr)
        sys.exit(1)
        
    try:
        df = pd.read_csv(file_path, sep=separator)
    except Exception as e:
        print(f"Error reading file: {e}", file=sys.stderr)
        sys.exit(1)

    if log2fc_col not in df.columns:
        print(f"Error: Column '{log2fc_col}' not found in the input file.", file=sys.stderr)
        print(f"Available columns are: {list(df.columns)}", file=sys.stderr)
        sys.exit(1)

    # --- 2. The Filtering Logic ---
    initial_gene_count = len(df)
    
    # The filtering "equation": keep rows where the log2FC is less than or equal to the threshold.
    # This removes genes with high positive fold changes, which are presumed to be contaminants.
    filtered_df = df[df[log2fc_col] <= threshold].copy()
    
    final_gene_count = len(filtered_df)
    removed_count = initial_gene_count - final_gene_count

    # --- 3. Report Summary and Output Data ---
    # Print the summary of the filtering operation to the standard error stream.
    # This prevents the summary from being mixed with the output data.
    print("--- Filtering Summary ---", file=sys.stderr)
    print(f"Filtering Rule Applied: Gene kept if '{log2fc_col}' <= {threshold}", file=sys.stderr)
    print(f"Initial number of genes: {initial_gene_count}", file=sys.stderr)
    print(f"Number of contaminant genes removed: {removed_count}", file=sys.stderr)
    print(f"Final number of genes: {final_gene_count}", file=sys.stderr)
    print("-------------------------", file=sys.stderr)
    print("Cleaned data is being printed to standard output...", file=sys.stderr)

    # Print the cleaned dataframe to standard output as a CSV.
    # This allows for redirection into a new file in a pipeline.
    filtered_df.to_csv(sys.stdout, index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Filter contaminant genes from a differential expression results file based on a log2 Fold Change (log2FC) threshold. This script assumes contaminants cause a large positive log2FC.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "input_file",
        type=str,
        help="Path to the input differential expression file (must be .csv or .tsv)."
    )
    parser.add_argument(
        "--log2fc_col",
        type=str,
        default="log2FoldChange",
        help="Name of the column containing log2 fold change values.\n(default: 'log2FoldChange')"
    )
    parser.add_argument(
        "--threshold",
        type=float,
        required=True,
        help="The log2FC threshold. Genes with a value in 'log2fc_col' GREATER than this will be removed.\n(e.g., 2.0 or 3.0)"
    )

    # Example of how to run from the command line:
    # python this_script_name.py /path/to/your/de_results.csv --threshold 2.5 > cleaned_de_results.csv
    #
    # This would read 'de_results.csv', remove all genes where 'log2FoldChange' > 2.5,
    # and save the new, cleaned dataset to 'cleaned_de_results.csv'.

    args = parser.parse_args()
    filter_contaminants_by_logfc(args.input_file, args.log2fc_col, args.threshold)