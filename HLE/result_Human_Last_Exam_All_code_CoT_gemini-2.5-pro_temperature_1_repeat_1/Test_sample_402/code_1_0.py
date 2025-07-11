import argparse
import csv
import sys

def filter_contaminants_by_logfc(input_file, logfc_col, threshold):
    """
    Filters a differential expression file to remove potential contaminants.

    Contaminants are identified as genes with a log2 fold change below a
    specified negative threshold. This script assumes the comparison was
    Test (CAR-IL15) vs. Control (CAR-only), and contamination is in the
    Control group.

    The filtered data is printed to standard output.
    """
    # Use csv.Sniffer to automatically detect if the file is comma or tab-separated
    try:
        with open(input_file, 'r', newline='') as f:
            dialect = csv.Sniffer().sniff(f.read(2048))
            f.seek(0)
            reader = csv.reader(f, dialect)
    except (csv.Error, FileNotFoundError) as e:
        print(f"Error reading input file: {e}", file=sys.stderr)
        sys.exit(1)

    # Read the header and find the column index for the logFC values
    header = next(reader)
    try:
        logfc_idx = header.index(logfc_col)
    except ValueError:
        print(
            f"Error: Column '{logfc_col}' not found in the input file header.",
            file=sys.stderr
        )
        print(f"Available columns are: {header}", file=sys.stderr)
        sys.exit(1)

    # Use a CSV writer for standard output, maintaining the original file format
    writer = csv.writer(sys.stdout, dialect)
    writer.writerow(header)

    # To satisfy the request "output each number in the final equation",
    # we print the exact filtering rule being applied to the error stream.
    print(f"Filtering Equation: Keep gene if value in column '{logfc_col}' >= {threshold}", file=sys.stderr)
    
    # Process the data rows
    initial_count = 0
    final_count = 0
    for row in reader:
        initial_count += 1
        try:
            # Skip rows with non-numeric or empty logFC values by keeping them in the output
            logfc_val_str = row[logfc_idx]
            if logfc_val_str in ('', 'NA', 'null'):
                writer.writerow(row)
                final_count += 1
                continue
                
            logfc_val = float(logfc_val_str)
            
            # This is the filtering condition based on the equation
            if logfc_val >= threshold:
                writer.writerow(row)
                final_count += 1

        except (ValueError, IndexError):
            print(f"Warning: Could not parse row: {row}. Skipping.", file=sys.stderr)
            continue

    # Print a final summary of the operation to the error stream
    removed_count = initial_count - final_count
    print(f"\nFiltering complete.", file=sys.stderr)
    print(f"Retained {final_count} of {initial_count} total genes.", file=sys.stderr)
    print(f"Removed {removed_count} potential contaminant genes.", file=sys.stderr)

def main():
    """Main function to parse command-line arguments and run the filter."""
    parser = argparse.ArgumentParser(
        description="Filter differential expression results to remove "
                    "contaminants based on a negative log2 fold change threshold.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Path to the input differential expression file (CSV or TSV format)."
    )
    parser.add_argument(
        "-c", "--column",
        required=True,
        help="Name of the log2 fold change column (e.g., 'log2FoldChange')."
    )
    parser.add_argument(
        "-t", "--threshold",
        type=float,
        default=-5.0,
        help="Log2 fold change threshold. Genes with a value BELOW this "
             "will be removed.\nDefault: -5.0"
    )
    args = parser.parse_args()

    filter_contaminants_by_logfc(args.input, args.column, args.threshold)

if __name__ == "__main__":
    main()