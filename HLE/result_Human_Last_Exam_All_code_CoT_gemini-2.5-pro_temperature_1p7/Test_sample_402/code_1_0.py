import pandas as pd
import argparse

def filter_contaminants(input_file, output_file, lfc_threshold, padj_threshold):
    """
    Filters a differential expression results file to remove contaminating genes.

    Contaminating genes are defined as those with a log2FoldChange below
    a given negative threshold and a significant adjusted p-value, assuming a
    comparison of (Test group) vs (Contaminated control group).

    Args:
        input_file (str): Path to the input CSV/TSV file.
        output_file (str): Path to save the filtered output CSV file.
        lfc_threshold (float): Negative log2 fold change threshold. Contaminants are below this.
        padj_threshold (float): Adjusted p-value threshold for significance.
    """
    # --- 1. Load the data ---
    try:
        # Auto-detect separator (CSV or TSV)
        sep = ',' if input_file.lower().endswith('.csv') else '\t'
        df = pd.read_csv(input_file, sep=sep)
        # Ensure required columns exist
        required_cols = ['log2FoldChange', 'padj']
        if not all(col in df.columns for col in required_cols):
            print(f"Error: Input file must contain the columns 'log2FoldChange' and 'padj'.")
            return
        print(f"Successfully loaded {input_file} containing {len(df)} total genes.")
    except FileNotFoundError:
        print(f"Error: Input file not found at {input_file}")
        return
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")
        return

    # --- 2. Define and apply the filter ---
    # Per your request, we will print the numbers used in the filtering equation.
    print("\n--- Filtering Strategy ---")
    print(f"Identifying contaminant genes using the rule:")
    print(f"log2FoldChange < {lfc_threshold} AND padj < {padj_threshold}")
    
    # Define the boolean condition that identifies a contaminant gene
    is_contaminant = (df['log2FoldChange'] < lfc_threshold) & (df['padj'] < padj_threshold)

    # Invert the condition to keep all genes that are NOT contaminants
    filtered_df = df[~is_contaminant].copy()
    contaminant_df = df[is_contaminant]

    # --- 3. Report Results ---
    num_total = len(df)
    num_removed = len(contaminant_df)
    num_remaining = len(filtered_df)

    print("\n--- Summary ---")
    print(f"Total genes before filtering: {num_total}")
    print(f"Number of contaminating genes removed: {num_removed}")
    print(f"Number of genes remaining after filtering: {num_remaining}")

    if num_removed > 0 and len(df.columns) > 0:
        print("\nTop 5 potential contaminating genes removed (sorted by lowest log2FoldChange):")
        gene_col_name = df.columns[0]
        print(contaminant_df.nsmallest(5, 'log2FoldChange')[[gene_col_name, 'log2FoldChange', 'padj']])

    # --- 4. Save the filtered data ---
    try:
        filtered_df.to_csv(output_file, index=False)
        print(f"\nFiltered data successfully saved to: {output_file}")
    except Exception as e:
        print(f"An error occurred while saving the file: {e}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Filter contaminating genes from RNA-seq differential expression results. "
                    "Assumes contamination is in the control group, resulting in a large negative log2FC.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '-i', '--input',
        type=str,
        required=True,
        help="Path to the input differential expression file (e.g., from DESeq2). Must be CSV or TSV.\n"
             "Expected columns: 'log2FoldChange', 'padj'."
    )
    parser.add_argument(
        '-o', '--output',
        type=str,
        required=True,
        help="Path for the output filtered CSV file."
    )
    parser.add_argument(
        '--lfc',
        type=float,
        default=-5.0,
        help="Log2 Fold Change threshold. Genes with LFC *below* this value will be flagged.\n"
             "A value of -5 implies the gene is 32-fold more abundant in the control.\n"
             "Default: -5.0"
    )
    parser.add_argument(
        '--padj',
        type=float,
        default=0.05,
        help="Adjusted p-value threshold. Genes must be below this to be considered significant.\n"
             "Default: 0.05"
    )

    args = parser.parse_args()

    # The LFC threshold must be negative for this logic to work
    if args.lfc >= 0:
        print("Error: The --lfc threshold must be a negative value to identify genes downregulated in the test group.")
    else:
        filter_contaminants(args.input, args.output, args.lfc, args.padj)