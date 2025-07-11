import sys
import argparse
import csv

def filter_contaminants_by_log2fc(input_file, output_file, log2fc_threshold, gene_col, log2fc_col, delimiter):
    """
    Filters a differential expression file to remove contaminants based on a log2FC threshold.

    Contaminants are assumed to be genes with a strong negative log2FC, indicating
    high expression in the control group but not the treatment group.
    """
    try:
        # Open the input file for reading and the output file for writing
        with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
            # Use the provided delimiter to create reader and writer objects
            reader = csv.reader(infile, delimiter=delimiter)
            writer = csv.writer(outfile, delimiter=delimiter)

            # Read the header and write it to the output file
            header = next(reader)
            writer.writerow(header)

            # Find the column indices for the gene name and log2FC
            try:
                gene_idx = header.index(gene_col)
                log2fc_idx = header.index(log2fc_col)
            except ValueError as e:
                print(f"Error: Column '{e}' not found in the header of {input_file}.", file=sys.stderr)
                print(f"Available columns are: {header}", file=sys.stderr)
                sys.exit(1)

            print(f"Filtering genes with {log2fc_col} <= {log2fc_threshold}", file=sys.stderr)
            filtered_count = 0
            
            # Process each gene row
            for row in reader:
                try:
                    # Skip empty rows
                    if not row:
                        continue

                    log2fc_value = float(row[log2fc_idx])

                    # The filtering condition
                    if log2fc_value > log2fc_threshold:
                        # If the gene passes the filter, write its full row to the output
                        writer.writerow(row)
                    else:
                        # If the gene is filtered, print a message to stderr and increment the counter
                        gene_name = row[gene_idx]
                        print(f"Filtered contaminant: {gene_name} (log2FC: {log2fc_value:.2f})", file=sys.stderr)
                        filtered_count += 1
                        
                except (ValueError, IndexError):
                    # Handle cases with malformed rows or non-numeric log2FC values
                    print(f"Warning: Could not process row: {row}", file=sys.stderr)
                    continue

            print(f"\nFiltering complete. Removed {filtered_count} potential contaminant genes.", file=sys.stderr)
            print(f"Filtered data saved to: {output_file}", file=sys.stderr)

    except FileNotFoundError:
        print(f"Error: The file '{input_file}' was not found.", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Filter potential contaminants from a differential expression file based on a negative log2 fold-change threshold. This is useful when the control group has specific contaminants not present in the treatment group.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        'input_file',
        help='Path to the input differential expression file (e.g., from DESeq2).'
    )
    parser.add_argument(
        'output_file',
        help='Path for the new filtered output file.'
    )
    parser.add_argument(
        '--threshold',
        type=float,
        default=-5.0,
        help='Log2 fold-change threshold. Genes with a value at or below this will be removed.\n'
             'A value of -5 means genes where expression in the control is >32-fold higher than\n'
             'in the treatment group will be filtered out. (default: -5.0)'
    )
    parser.add_argument(
        '--gene_col',
        type=str,
        default='gene_name',
        help="Name of the column containing gene identifiers. (default: 'gene_name')"
    )
    parser.add_argument(
        '--log2fc_col',
        type=str,
        default='log2FoldChange',
        help="Name of the log2 fold-change column. (default: 'log2FoldChange')"
    )
    parser.add_argument(
        '--delimiter',
        type=str,
        default='\t',
        help="Delimiter used in the input file. Use ',' for CSV and '\\t' for TSV. (default: '\\t')"
    )

    args = parser.parse_args()

    # Convert '\\t' from command line to an actual tab character
    delimiter_char = args.delimiter.replace('\\t', '\t')

    filter_contaminants_by_log2fc(
        input_file=args.input_file,
        output_file=args.output_file,
        log2fc_threshold=args.threshold,
        gene_col=args.gene_col,
        log2fc_col=args.log2fc_col,
        delimiter=delimiter_char
    )