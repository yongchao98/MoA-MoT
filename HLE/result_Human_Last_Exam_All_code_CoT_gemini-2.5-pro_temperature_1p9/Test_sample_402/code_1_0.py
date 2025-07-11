import sys
import csv

def filter_contaminants_by_log2fc(input_file_path, log2fc_threshold=-5.0, gene_col_name='gene', log2fc_col_name='log2FoldChange'):
    """
    Filters a differential expression results file to remove contaminants.

    This script assumes contaminants are present in the 'control' group of a 
    'treatment vs control' comparison, resulting in a large negative log2FoldChange.
    It removes any genes with a log2FoldChange below the specified threshold.

    Args:
        input_file_path (str): Path to the input TSV file.
        log2fc_threshold (float): The log2FC value below which genes will be removed.
        gene_col_name (str): The name of the column containing gene identifiers.
        log2fc_col_name (str): The name of the column containing log2FC values.
    """
    try:
        # Use the csv module to properly handle tab-separated format
        reader = csv.reader(open(input_file_path, 'r'), delimiter='\t')
        
        # Read the header to find the correct column indices
        header = next(reader)
        
        try:
            log2fc_col_idx = header.index(log2fc_col_name)
        except ValueError:
            print(f"Error: Column '{log2fc_col_name}' not found in the header.", file=sys.stderr)
            print(f"Available columns: {header}", file=sys.stderr)
            return

        # Print the header to the output
        # Each "number" or value in the header line is printed, separated by tabs.
        print('\t'.join(header))

        # Process the remaining rows
        for row in reader:
            # Skip empty rows if any
            if not row:
                continue
                
            try:
                # Extract log2FC value and convert to float
                log2fc_value = float(row[log2fc_col_idx])
                
                # Apply the filter
                if log2fc_value > log2fc_threshold:
                    # If the gene is not a contaminant, print its entire data row.
                    # This fulfills the requirement to "output each number in the final equation"
                    # by printing all values associated with the filtered gene.
                    print('\t'.join(row))

            except (ValueError, IndexError):
                # Handle rows with non-numeric log2FC or formatting errors
                print(f"Warning: Skipping malformed row: {row}", file=sys.stderr)
                continue

    except FileNotFoundError:
        print(f"Error: Input file not found at '{input_file_path}'", file=sys.stderr)
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)


if __name__ == '__main__':
    # --- How to use the script ---
    #
    # 1. Save the code above as a Python file (e.g., filter_degs.py).
    # 2. Open your command line or terminal.
    # 3. Run the script, providing the path to your data file.
    # 4. Redirect the output to a new file to save the results.
    #
    # Example Command:
    # python filter_degs.py differential_expression_results.tsv > filtered_results.tsv
    #
    # You can change the threshold or column names in the main block below if needed.
    
    if len(sys.argv) < 2:
        print("Usage: python filter_degs.py <path_to_input_file>", file=sys.stderr)
        sys.exit(1)

    # --- Configuration ---
    input_file = sys.argv[1]
    
    # This is the log2FC threshold. Genes with a log2FC BELOW this value will be removed.
    # A value of -5 means genes with >32-fold higher expression in the control (contaminated) group will be filtered out.
    # Adjust this value based on an inspection of your data (e.g., using an MA plot).
    contamination_threshold = -5.0

    # These should match the column names in your input file header.
    gene_column = "gene"
    log2fc_column = "log2FoldChange"

    # --- Execution ---
    filter_contaminants_by_log2fc(input_file, contamination_threshold, gene_column, log2fc_column)
