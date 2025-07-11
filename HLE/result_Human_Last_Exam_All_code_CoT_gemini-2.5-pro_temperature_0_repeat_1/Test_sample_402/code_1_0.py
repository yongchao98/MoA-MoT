import sys
import csv
import argparse

def filter_by_log2fc(input_file, output_file, gene_col, log2fc_col, threshold):
    """
    Filters a differential expression file to remove genes based on a log2FC threshold.

    This is useful for removing contaminant signatures that appear as strongly
    downregulated genes.
    """
    # Print the filtering parameters being used to standard error for user information.
    # This fulfills the requirement to output the numbers in the filtering equation.
    print(f"Filtering strategy: Keep genes where '{log2fc_col}' > {threshold}", file=sys.stderr)
    
    try:
        reader = csv.reader(input_file, delimiter='\t')
        writer = csv.writer(output_file, delimiter='\t', lineterminator='\n')

        # Find column indices from header
        header = next(reader)
        writer.writerow(header)

        try:
            gene_idx = header.index(gene_col)
            log2fc_idx = header.index(log2fc_col)
        except ValueError as e:
            print(f"Error: Column '{e.args[0]}' not found in the input file header.", file=sys.stderr)
            print(f"Available columns are: {header}", file=sys.stderr)
            sys.exit(1)

        # Process each row
        kept_count = 0
        removed_count = 0
        for row in reader:
            try:
                # Extract log2FC value and convert to float
                log2fc_value = float(row[log2fc_idx])

                # The filtering equation: check if the value is above the threshold
                if log2fc_value > threshold:
                    writer.writerow(row)
                    kept_count += 1
                else:
                    removed_count += 1
            except (ValueError, IndexError):
                # Skip rows with non-numeric or missing log2FC values
                print(f"Warning: Skipping malformed row: {row}", file=sys.stderr)
                continue
        
        print(f"\nFiltering complete.", file=sys.stderr)
        print(f"Kept {kept_count} genes.", file=sys.stderr)
        print(f"Removed {removed_count} genes based on the threshold.", file=sys.stderr)

    except FileNotFoundError:
        print(f"Error: Input file not found at {input_file.name}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Filter differential expression results based on a log2 fold change threshold. "
                    "Reads from a tab-delimited file and prints the filtered output to standard output.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "input_file",
        type=argparse.FileType('r'),
        help="Path to the input differential expression file (tab-separated)."
    )
    parser.add_argument(
        "--gene_col",
        type=str,
        default="gene",
        help="Name of the column containing gene identifiers (default: 'gene')."
    )
    parser.add_argument(
        "--log2fc_col",
        type=str,
        default="log2FoldChange",
        help="Name of the column containing log2 fold change values (default: 'log2FoldChange')."
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=-3.0,
        help="The log2FC threshold. Genes with a log2FC *below* this value will be removed (default: -3.0)."
    )

    args = parser.parse_args()

    # The script writes the filtered data to standard output (sys.stdout)
    # and informational messages to standard error (sys.stderr).
    filter_by_log2fc(
        input_file=args.input_file,
        output_file=sys.stdout,
        gene_col=args.gene_col,
        log2fc_col=args.log2fc_col,
        threshold=args.threshold
    )
