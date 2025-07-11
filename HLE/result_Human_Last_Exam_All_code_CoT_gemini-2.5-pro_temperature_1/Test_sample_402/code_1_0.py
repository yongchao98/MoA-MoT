import sys
import csv
import os

def filter_contaminating_genes():
    """
    Filters genes from a differential expression results file based on a log2FoldChange threshold.

    This script assumes that contaminating genes are highly expressed in the first condition
    of the comparison (e.g., 'CAR_only' in 'CAR_only vs. CAR_IL15'), resulting in a large
    positive log2FoldChange.
    """
    # --- Parameters ---
    # The input file should be a tab-separated file from your analysis.
    input_filename = 'differential_expression_results.tsv'
    
    # The column name for gene identifiers.
    gene_id_col = 'gene_id'
    
    # The column from which to read the log2 fold change values.
    # This assumes the comparison was of the form: contaminated_group vs. clean_group
    log2fc_col = 'log2FoldChange'
    
    # Set the filtering threshold. Genes with a log2FoldChange ABOVE this
    # value will be removed. A value of 2.0 means we filter genes that are
    # over 4-fold more expressed in the contaminated samples.
    log2fc_threshold = 2.0

    # --- Pre-run Check ---
    # For demonstration, create a sample input file if it doesn't exist.
    if not os.path.exists(input_filename):
        print(
            f"INFO: Sample file '{input_filename}' not found. Creating a demo file.",
            file=sys.stderr
        )
        with open(input_filename, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['gene_id', 'baseMean', 'log2FoldChange', 'pvalue', 'padj'])
            writer.writerow(['GENE_A', '150.5', '-1.58', '0.041', '0.09'])
            writer.writerow(['GENE_B', '275.1', '1.85', '0.012', '0.03'])
            writer.writerow(['CONTAMINANT_1', '890.2', '4.50', '0.0001', '0.0004'])
            writer.writerow(['GENE_C', '98.0', '0.50', '0.650', '0.78'])
            writer.writerow(['CONTAMINANT_2', '650.7', '3.15', '0.0005', '0.0015'])
            writer.writerow(['GENE_D', '432.9', '1.99', '0.009', '0.025'])

    # --- Main Logic ---
    try:
        with open(input_filename, 'r', newline='') as infile:
            # Use DictReader to easily access columns by name.
            reader = csv.DictReader(infile, delimiter='\t')
            
            # The writer will print to standard output (the console).
            # This allows you to redirect the output to a new file.
            writer = csv.DictWriter(sys.stdout, fieldnames=reader.fieldnames, delimiter='\t')
            
            # Print the header row to the output.
            writer.writeheader()
            
            # Iterate over each gene in the input file.
            for row in reader:
                try:
                    log2fc_value = float(row[log2fc_col])
                    
                    # This is the filtering equation.
                    # We KEEP the gene if its log2FC is not greater than the threshold.
                    if log2fc_value <= log2fc_threshold:
                        writer.writerow(row)
                    else:
                        # For filtered genes, print the equation and values to stderr.
                        gene_id = row[gene_id_col]
                        print(
                            f"FILTERED: Gene '{gene_id}' was removed because its expression change "
                            f"({log2fc_value:.2f}) exceeded the threshold ({log2fc_threshold}). "
                            f"Equation: {log2fc_value:.2f} > {log2fc_threshold}",
                            file=sys.stderr
                        )

                except (ValueError, TypeError):
                    print(f"WARNING: Could not parse log2FC value for row: {row}. Skipping.", file=sys.stderr)
                except KeyError:
                    print(f"ERROR: Required column '{log2fc_col}' or '{gene_id_col}' not found in the file.", file=sys.stderr)
                    sys.exit(1)
                    
    except FileNotFoundError:
        print(f"ERROR: Input file '{input_filename}' could not be read.", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    filter_contaminating_genes()