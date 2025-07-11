import pandas as pd
import argparse
import sys

def filter_contaminants_by_logfc(input_file, output_file, logfc_column, threshold):
    """
    Filters out potential contaminant genes based on a negative log2 fold change threshold.

    Args:
        input_file (str): Path to the input differential expression results file (CSV/TSV).
        output_file (str): Path to save the filtered output file.
        logfc_column (str): The name of the column containing log2 fold change values.
        threshold (float): The negative log2FC threshold. Genes with a log2FC below
                           this value will be removed.
    """
    # --- 1. Parameter Validation ---
    if threshold > 0:
        print(f"Error: The threshold must be a negative value to filter for contaminants in the control group. Provided: {threshold}", file=sys.stderr)
        sys.exit(1)

    # --- 2. Read Input Data ---
    try:
        # Use a flexible separator to handle both CSV and TSV files
        df = pd.read_csv(input_file, sep=None, engine='python')
        if logfc_column not in df.columns:
            print(f"Error: The specified logFC column '{logfc_column}' was not found in the input file.", file=sys.stderr)
            print(f"Available columns are: {list(df.columns)}", file=sys.stderr)
            sys.exit(1)
    except FileNotFoundError:
        print(f"Error: Input file not found at '{input_file}'", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred while reading the file: {e}", file=sys.stderr)
        sys.exit(1)

    initial_gene_count = len(df)
    print(f"Successfully loaded {initial_gene_count} genes from '{input_file}'.")

    # --- 3. Apply Filtering Logic ---
    # The equation for filtering is: Keep genes where log2FoldChange >= threshold
    # The final output will print this equation with the numbers used.
    print("-" * 50)
    print("Applying filtering based on the following equation:")
    print(f"Keep gene if: value_in_column('{logfc_column}') >= {threshold}")
    print("-" * 50)
    
    # Create a boolean mask for genes to keep
    genes_to_keep_mask = df[logfc_column] >= threshold
    
    # Separate the data into kept and removed dataframes
    filtered_df = df[genes_to_keep_mask]
    removed_df = df[~genes_to_keep_mask]
    
    genes_removed_count = len(removed_df)
    final_gene_count = len(filtered_df)

    # --- 4. Report Results ---
    print(f"Identified and removed {genes_removed_count} potential contaminant genes.")
    if genes_removed_count > 0:
        print("\n--- Removed Genes (Top 5) ---")
        print(removed_df.head())
        print("-" * 30)

    print(f"\nTotal genes remaining: {final_gene_count}")

    # --- 5. Save Filtered Data ---
    try:
        filtered_df.to_csv(output_file, index=False)
        print(f"Filtered data successfully saved to '{output_file}'.")
    except Exception as e:
        print(f"An error occurred while saving the file: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    # Set up the argument parser for command-line execution
    parser = argparse.ArgumentParser(
        description="Filter contaminant genes from RNA-seq differential expression results based on a log2 Fold Change threshold.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '-i', '--input',
        required=True,
        help="Path to the input differential expression file (e.g., results.csv)."
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help="Path for the output filtered file (e.g., filtered_results.csv)."
    )
    parser.add_argument(
        '-c', '--column',
        default='log2FoldChange',
        help="Name of the log2 Fold Change column in the input file (default: 'log2FoldChange')."
    )
    parser.add_argument(
        '-t', '--threshold',
        type=float,
        default=-5.0,
        help="Negative log2FC threshold. Genes with a log2FC strictly BELOW this value will be removed (default: -5.0)."
    )

    args = parser.parse_args()

    # To run this script, save it as `filter_contaminants.py` and execute from your terminal:
    #
    # Example Usage:
    # python filter_contaminants.py --input your_de_results.csv --output your_filtered_results.csv --threshold -5
    #
    # To test, you can create a dummy CSV file named 'your_de_results.csv' with the following content:
    #
    # GeneName,log2FoldChange,pvalue,padj
    # GENE_A,3.5,0.001,0.005
    # GENE_B,-1.2,0.2,0.4
    # TUMOR_GENE_1,-8.1,0.00001,0.00009
    # GENE_C,0.5,0.6,0.8
    # TUMOR_GENE_2,-6.5,0.00004,0.0002
    #

    filter_contaminants_by_logfc(
        input_file=args.input,
        output_file=args.output,
        logfc_column=args.column,
        threshold=args.threshold
    )