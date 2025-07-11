import sys
import argparse
import csv

def filter_contaminants():
    """
    Filters a differential expression results file to remove potential contaminants
    based on a log2 Fold Change threshold.
    """
    parser = argparse.ArgumentParser(
        description="Filter out contaminating genes from a differential expression results file. "
                    "This script assumes contaminants are highly expressed in the control group, "
                    "resulting in a large negative log2FoldChange. It prints the filtered data to standard output.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "input_file",
        help="Path to the input differential expression file (e.g., 'DE_results.csv')."
    )
    parser.add_argument(
        "--lfc_col",
        default="log2FoldChange",
        help="Name of the log2 Fold Change column in your file (default: 'log2FoldChange')."
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=-2.0,
        help="Log2FC threshold. Genes with a value BELOW this will be removed (default: -2.0)."
    )

    args = parser.parse_args()

    try:
        with open(args.input_file, 'r', newline='') as infile:
            # Sniff the file to automatically detect the delimiter (comma or tab)
            try:
                dialect = csv.Sniffer().sniff(infile.read(2048))
            except csv.Error:
                # Default to tab-separated if sniffing fails
                dialect = csv.excel_tab
            infile.seek(0)
            reader = csv.reader(infile, dialect)
            
            header = next(reader)
            
            # Find the index of the log2FoldChange column
            try:
                lfc_idx = header.index(args.lfc_col)
            except ValueError:
                print(f"Error: Column '{args.lfc_col}' not found in the file header.", file=sys.stderr)
                print(f"Available columns are: {', '.join(header)}", file=sys.stderr)
                sys.exit(1)

            # Print the header for the output file
            print(dialect.delimiter.join(header))
            
            total_rows = 0
            kept_rows = 0
            
            # The filtering equation is: gene_log2FC >= threshold
            for row in reader:
                total_rows += 1
                try:
                    # Get the log2FC value for the current gene
                    log2fc_val = float(row[lfc_idx])
                    
                    # Apply the filtering condition
                    if log2fc_val >= args.threshold:
                        # If the condition is met, print the entire row (all its numbers)
                        print(dialect.delimiter.join(row))
                        kept_rows += 1
                except (ValueError, IndexError):
                    # This handles cases of non-numeric or missing LFC values
                    print(f"Warning: Could not parse log2FC value in row: {row}", file=sys.stderr)
                    continue
        
        # Print a summary of the filtering process to the standard error stream
        print(f"\n--- Filtering Summary ---", file=sys.stderr)
        print(f"Filtering complete for file: {args.input_file}", file=sys.stderr)
        print(f"Filtering condition: Genes kept if '{args.lfc_col}' >= {args.threshold}", file=sys.stderr)
        print(f"Total genes processed: {total_rows}", file=sys.stderr)
        print(f"Genes kept: {kept_rows}", file=sys.stderr)
        print(f"Genes filtered out (potential contaminants): {total_rows - kept_rows}", file=sys.stderr)

    except FileNotFoundError:
        print(f"Error: Input file not found at '{args.input_file}'", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    filter_contaminants()