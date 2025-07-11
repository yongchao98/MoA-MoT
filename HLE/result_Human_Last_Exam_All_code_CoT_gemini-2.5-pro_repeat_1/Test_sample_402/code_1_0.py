import pandas as pd
import sys

# --- Parameters to Modify ---

# 1. Path to your differential expression results file.
#    The script assumes a CSV or TSV format.
INPUT_FILE = 'de_results.csv'

# 2. Names of the relevant columns in your file (case-sensitive).
GENE_COLUMN = 'gene_name'        # Column with gene identifiers
LOG2FC_COLUMN = 'log2FoldChange' # Column with log2 fold change values
PADJ_COLUMN = 'padj'             # Column with adjusted p-value

# 3. Filtering thresholds.
#    The script will REMOVE genes that satisfy both conditions:
#    (log2FoldChange < LOG2FC_THRESHOLD) AND (padj < PADJ_THRESHOLD)
LOG2FC_THRESHOLD = -2.0
PADJ_THRESHOLD = 0.05

# --- Script Logic (No modifications needed below this line) ---

def filter_contaminants():
    """
    Reads a differential expression file, filters out contaminating genes
    based on log2FC and padj, and prints the result.
    """
    try:
        # Attempt to read the file, trying comma and then tab separators.
        try:
            df = pd.read_csv(INPUT_FILE)
        except (pd.errors.ParserError, UnicodeDecodeError):
            df = pd.read_csv(INPUT_FILE, sep='\t')
    except FileNotFoundError:
        print(f"Error: The file '{INPUT_FILE}' was not found.", file=sys.stderr)
        print("Please update the 'INPUT_FILE' variable in the script with the correct path.", file=sys.stderr)
        return
    except Exception as e:
        print(f"An error occurred while reading the file: {e}", file=sys.stderr)
        return

    # Check if required columns exist
    required_cols = [LOG2FC_COLUMN, PADJ_COLUMN]
    for col in required_cols:
        if col not in df.columns:
            print(f"Error: Column '{col}' not found in the input file.", file=sys.stderr)
            print("Please check the column name variables (LOG2FC_COLUMN, PADJ_COLUMN) in the script.", file=sys.stderr)
            return

    # Define the filtering equation for identifying contaminants
    is_contaminant = (df[LOG2FC_COLUMN] < LOG2FC_THRESHOLD) & (df[PADJ_COLUMN] < PADJ_THRESHOLD)

    # Separate the dataframe into kept genes and removed contaminants
    original_gene_count = len(df)
    df_filtered = df[~is_contaminant].copy()
    contaminant_count = original_gene_count - len(df_filtered)

    # Print the filtering summary
    print("--- Filtering Summary ---")
    print(f"Filtering based on the following equation to identify and remove contaminants:")
    print(f"    {LOG2FC_COLUMN} < {LOG2FC_THRESHOLD}")
    print(f"    AND")
    print(f"    {PADJ_COLUMN} < {PADJ_THRESHOLD}")
    print("-" * 25)
    print(f"Original number of genes: {original_gene_count}")
    print(f"Number of contaminant genes removed: {contaminant_count}")
    print(f"Final number of genes after filtering: {len(df_filtered)}")
    print("-" * 25)
    print("\n--- Filtered Data ---")
    
    # Print the cleaned dataframe to the console.
    # To save this output to a file, run the script from your terminal like this:
    # python your_script_name.py > filtered_results.csv
    print(df_filtered.to_string())

if __name__ == '__main__':
    filter_contaminants()