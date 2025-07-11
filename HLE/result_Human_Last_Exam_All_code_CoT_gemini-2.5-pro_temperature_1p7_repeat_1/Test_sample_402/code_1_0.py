import argparse
import pandas as pd
import sys

def filter_de_results(input_file, output_file, lfc_col, threshold):
    """
    Filters a differential expression results file to remove genes
    with a log2FoldChange below a specified negative threshold.

    This is designed to remove contaminants assumed to be present in the
    'control' condition of a 'treatment vs control' comparison.
    """
    print(f"--- Contaminant Filtering Script ---")

    # --- Step 1: Read the input file ---
    try:
        # Detect separator based on file extension
        if input_file.endswith('.tsv'):
            sep = '\t'
        else:
            sep = ','  # Default to comma for .csv or other extensions
        
        df = pd.read_csv(input_file, sep=sep)
        print(f"Successfully read {len(df)} total genes from '{input_file}'.")
    except FileNotFoundError:
        print(f"Error: The input file '{input_file}' was not found.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: Failed to read or parse the input file. Please ensure it is a valid CSV or TSV file.", file=sys.stderr)
        print(f"Details: {e}", file=sys.stderr)
        sys.exit(1)

    # --- Step 2: Validate that the log2FC column exists ---
    if lfc_col not in df.columns:
        print(f"Error: The specified log2FoldChange column '{lfc_col}' was not found in the file.", file=sys.stderr)
        print(f"Available columns are: {list(df.columns)}", file=sys.stderr)
        sys.exit(1)

    # --- Step 3: Apply the filter ---
    # The filtering "equation" is: log2FoldChange >= threshold
    # We print the components of this equation for clarity.
    print("\nFiltering Equation:")
    print(f"  Gene attribute: '{lfc_col}'")
    print(f"  Operator: >=")
    print(f"  Threshold value: {threshold}")
    
    # Perform the filtering
    filtered_df = df[df[lfc_col] >= threshold].copy()
    
    num_original = len(df)
    num_filtered = len(filtered_df)
    num_removed = num_original - num_filtered
    
    print(f"\nFiltering complete.")
    print(f"  {num_removed} genes were removed based on the threshold.")
    print(f"  {num_filtered} genes remain in the filtered dataset.")

    # --- Step 4: Write the filtered data to the output file ---
    try:
        # Use the same separator for writing as was used for reading
        filtered_df.to_csv(output_file, sep=sep, index=False)
        print(f"\nSuccessfully wrote filtered results to '{output_file}'.")
    except Exception as e:
        print(f"Error: Could not write to output file '{output_file}'.", file=sys.stderr)
        print(f"Details: {e}", file=sys.stderr)
        sys.exit(1)
        
    print("--- Script Finished ---")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Filter differential expression (DE) results to remove potential contaminants. "
                    "This script removes genes with a log2FoldChange below a specified negative threshold.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    parser.add_argument(
        '-i', '--input', 
        required=True, 
        help="Path to the input DE results file (CSV or TSV format)."
    )
    
    parser.add_argument(
        '-o', '--output', 
        required=True, 
        help="Path for the filtered output file."
    )
    
    parser.add_argument(
        '--lfc_col', 
        default='log2FoldChange', 
        help="Name of the log2FoldChange column in the input file.\n(Default: 'log2FoldChange')"
    )
    
    parser.add_argument(
        '-t', '--threshold', 
        type=float, 
        default=-5.0,
        help="The log2FoldChange threshold. Genes with a value *below* this will be removed.\n(Default: -5.0, removes genes with LFC < -5)"
    )
    
    args = parser.parse_args()
    
    # Ensure threshold is negative, as that's the logic for this problem
    if args.threshold > 0:
        print(f"Warning: The provided threshold {args.threshold} is positive. This script is designed to remove contaminants "
              f"with large *negative* log2FoldChanges. Proceeding, but ensure this is your intent.", file=sys.stderr)

    filter_de_results(
        input_file=args.input,
        output_file=args.output,
        lfc_col=args.lfc_col,
        threshold=args.threshold
    )
    
    # Example usage from your command line:
    # python this_script_name.py \
    #   -i /path/to/your/de_results.csv \
    #   -o /path/to/your/filtered_results.csv \
    #   -t -5.0 \
    #   --lfc_col "log2FoldChange"
