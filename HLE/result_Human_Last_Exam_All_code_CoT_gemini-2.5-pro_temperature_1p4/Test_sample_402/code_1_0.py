import pandas as pd
import sys

def filter_contaminant_genes():
    """
    Reads a differential expression file, filters out contaminant genes based on
    log2FoldChange and p-value, and prints the cleaned data to standard output.
    """
    # =========================================================================
    # ## USER CONFIGURATION ##
    # Please edit the following parameters to match your setup.
    # =========================================================================

    # 1. Path to your differential expression results file.
    # This file should contain the comparison of (CAR T + IL15) vs (CAR T only).
    input_filepath = "path/to/your/de_results.csv"

    # 2. Names of the relevant columns in your file.
    # Check your file header for the correct names (e.g., 'logFC', 'adj.P.Val').
    log2fc_column = "log2FoldChange"
    padj_column = "padj"

    # 3. Filtering thresholds.
    # A gene is considered a contaminant if its log2FC AND padj are below these values.
    # A log2FC of -2.0 means a 4-fold lower expression in the IL15 group.
    log2fc_threshold = -2.0
    padj_threshold = 0.05
    # =========================================================================
    
    # Explain the filtering equation being used
    print("Filtering Strategy Report:", file=sys.stderr)
    print("----------------------------------------------------------------", file=sys.stderr)
    print("The final equation to identify and remove a contaminating gene is:", file=sys.stderr)
    print(f"({log2fc_column} < {log2fc_threshold}) AND ({padj_column} < {padj_threshold})", file=sys.stderr)
    print("----------------------------------------------------------------\n", file=sys.stderr)

    try:
        # Read the file, assuming the first column contains gene names/IDs
        sep = '\t' if input_filepath.lower().endswith('.tsv') else ','
        de_data = pd.read_csv(input_filepath, sep=sep, index_col=0)
    except FileNotFoundError:
        print(f"Error: The file '{input_filepath}' was not found.", file=sys.stderr)
        print("Please update the 'input_filepath' variable in the script.", file=sys.stderr)
        return
    except Exception as e:
        print(f"Error reading the data file: {e}", file=sys.stderr)
        return

    # Validate that the specified columns exist
    if log2fc_column not in de_data.columns or padj_column not in de_data.columns:
        print(f"Error: One or both specified columns ('{log2fc_column}', '{padj_column}') were not found.", file=sys.stderr)
        print(f"Available columns in your file are: {list(de_data.columns)}", file=sys.stderr)
        print("Please update the column name variables in the script.", file=sys.stderr)
        return

    print(f"Original gene count: {len(de_data)}", file=sys.stderr)

    # Ensure data types are numeric, converting non-numeric values to NaN (Not a Number)
    de_data[log2fc_column] = pd.to_numeric(de_data[log2fc_column], errors='coerce')
    de_data[padj_column] = pd.to_numeric(de_data[padj_column], errors='coerce')

    # Create the filtering condition for contaminant genes
    is_contaminant = (de_data[log2fc_column] < log2fc_threshold) & \
                     (de_data[padj_column] < padj_threshold)
                     
    # The '.fillna(False)' handles cases where p-value or log2FC might be NaN
    is_contaminant = is_contaminant.fillna(False)

    # Separate the clean data from the contaminants
    cleaned_data = de_data[~is_contaminant]
    
    genes_removed_count = is_contaminant.sum()
    print(f"Contaminating genes removed: {genes_removed_count}", file=sys.stderr)
    print(f"Final gene count after filtering: {len(cleaned_data)}\n", file=sys.stderr)
    print("Writing cleaned data to standard output...", file=sys.stderr)

    # Print the final cleaned dataframe to standard output in TSV format
    # To save this to a file, run the script like this:
    # python your_script_name.py > cleaned_results.tsv
    cleaned_data.to_csv(sys.stdout, sep='\t')


if __name__ == "__main__":
    filter_contaminant_genes()