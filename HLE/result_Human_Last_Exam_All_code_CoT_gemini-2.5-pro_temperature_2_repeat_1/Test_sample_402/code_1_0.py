import argparse
import pandas as pd
import sys
import os

def filter_contaminants_by_l2fc(input_file, l2fc_threshold, gene_col, l2fc_col):
    """
    Filters a differential expression results file to identify potential contaminant genes
    based on a log2 fold change threshold.

    The "equation" for filtering is: log2FoldChange > threshold.
    This script prints the genes that satisfy this condition.

    Args:
        input_file (str): Path to the CSV or TSV file from a DE analysis.
        l2fc_threshold (float): The log2 fold change threshold for filtering.
        gene_col (str): The name of the column containing gene identifiers.
        l2fc_col (str): The name of the column containing log2 fold change values.
    """
    # --- 1. Validate Input ---
    if not os.path.exists(input_file):
        print(f"Error: Input file not found at '{input_file}'", file=sys.stderr)
        sys.exit(1)

    # --- 2. Read Data ---
    try:
        # Auto-detect separator (CSV or TSV)
        sep = ',' if input_file.endswith('.csv') else '\t'
        de_results = pd.read_csv(input_file, sep=sep)
    except Exception as e:
        print(f"Error reading the input file: {e}", file=sys.stderr)
        sys.exit(1)

    # --- 3. Validate Columns ---
    if gene_col not in de_results.columns or l2fc_col not in de_results.columns:
        print(f"Error: Required columns '{gene_col}' or '{l2fc_col}' not found in the file.", file=sys.stderr)
        print(f"Available columns: {list(de_results.columns)}", file=sys.stderr)
        sys.exit(1)
        
    # --- 4. Filter for Contaminants ---
    # The filtering condition: find genes where the log2FoldChange is greater than the threshold.
    # These are genes highly enriched in the contaminated group vs. the clean group.
    contaminant_genes = de_results[de_results[l2fc_col] > l2fc_threshold].copy()
    
    # --- 5. Output Results ---
    # To satisfy the "output each number in the final equation" requirement,
    # we print the gene name, its L2FC value (the number), the comparison operator (>),
    # and the threshold (the other number).
    print(f"Filtering for genes where {l2fc_col} > {l2fc_threshold}")
    print("-" * 40)
    
    if contaminant_genes.empty:
        print("No contaminant genes found with the given threshold.")
    else:
        print("Potential Contaminant Genes to Remove:")
        # We print the gene name and the L2FC value that passed the filter.
        # This output can be redirected to a file for later use.
        output_df = contaminant_genes[[gene_col, l2fc_col]]
        # Print to stdout in a clean, machine-readable format (CSV)
        print(output_df.to_csv(index=False))

if __name__ == '__main__':
    # --- How to use this script from the command line ---
    #
    # 1. Save this code as a Python file, e.g., "filter_contaminants.py".
    #
    # 2. Create a dummy differential expression results file named "de_results.csv":
    #    For example, in your shell:
    #    cat > de_results.csv <<EOF
    # gene_symbol,baseMean,log2FoldChange,pvalue
    # GENE_A,1500,5.6,0.0001
    # GENE_B,250,1.2,0.04
    # CANCER_GENE_1,5000,8.2,0.00001
    # GENE_C,100,-2.5,0.001
    # CANCER_GENE_2,2100,6.5,0.00009
    # EOF
    #
    # 3. Run the script from your terminal:
    #    python filter_contaminants.py --input_file de_results.csv --l2fc_threshold 4.0
    #
    # 4. To save the list for the next step in your pipeline:
    #    python filter_contaminants.py --input_file de_results.csv --l2fc_threshold 4.0 > contaminant_genes.csv
    #
    parser = argparse.ArgumentParser(description="Filter contaminant genes from DE results based on L2FC.")
    parser.add_argument('--input_file', type=str, required=True, help="Path to the DE results file (CSV or TSV).")
    parser.add_argument('--l2fc_threshold', type=float, default=2.0, help="Log2 fold change threshold. Genes above this will be filtered.")
    parser.add_argument('--gene_col', type=str, default='gene_symbol', help="Name of the column with gene identifiers.")
    parser.add_argument('--l2fc_col', type=str, default='log2FoldChange', help="Name of the column with log2 fold change values.")
    
    args = parser.parse_args()
    
    filter_contaminants_by_l2fc(args.input_file, args.l2fc_threshold, args.gene_col, args.l2fc_col)
