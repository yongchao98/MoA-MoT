import argparse
import pandas as pd

def filter_contaminants(input_file, output_file, log2fc_threshold, padj_threshold):
    """
    Filters a differential expression results file to identify and extract
    contaminating genes based on log2 fold change and adjusted p-value.

    The filtering logic assumes that the 'contaminated' group is the numerator
    in the log2 fold change calculation (e.g., 'Contaminated_vs_Clean').
    Contaminating genes will therefore have a high positive log2FoldChange.

    Args:
        input_file (str): Path to the input CSV/TSV file from a DE analysis tool.
                          Expected columns: 'gene_id', 'log2FoldChange', 'padj'.
        output_file (str): Path to save the list of contaminating gene IDs.
        log2fc_threshold (float): The minimum log2 fold change to be considered a contaminant.
        padj_threshold (float): The maximum adjusted p-value to be considered significant.
    """
    try:
        # Detect separator (CSV or TSV) and read the data
        with open(input_file, 'r') as f:
            line = f.readline()
            if '\t' in line:
                sep = '\t'
            else:
                sep = ','
        
        de_results = pd.read_csv(input_file, sep=sep)
        
        # Ensure required columns exist
        required_cols = ['log2FoldChange', 'padj']
        # Try to find a gene identifier column. Common names are checked.
        gene_col = None
        for col_name in ['gene', 'gene_id', 'id', 'gene_name', de_results.columns[0]]:
            if col_name in de_results.columns:
                gene_col = col_name
                break
        
        if not all(col in de_results.columns for col in required_cols) or not gene_col:
            raise ValueError(f"Input file must contain a gene identifier column and the columns 'log2FoldChange' and 'padj'.")

        print(f"Successfully loaded {input_file} with {len(de_results)} total genes.")
        print("Using gene identifier column:", gene_col)
        
        # Define the filtering equation based on thresholds
        # Contaminants have high positive log2FC and low padj
        filter_condition = (de_results['log2FoldChange'] > log2fc_threshold) & (de_results['padj'] < padj_threshold)
        
        # Apply the filter
        contaminant_genes_df = de_results[filter_condition]
        
        # Extract the gene IDs
        contaminant_gene_list = contaminant_genes_df[gene_col]
        
        # Save the list of contaminant genes to the output file
        contaminant_gene_list.to_csv(output_file, index=False, header=False)
        
        print("\n--- Filtering Complete ---")
        print(f"Identified {len(contaminant_gene_list)} potential contaminating genes.")
        
        # Output the final filtering equation with the numbers used
        print("\nFinal filtering equation applied:")
        print(f"log2FoldChange > {log2fc_threshold}")
        print(f"padj < {padj_threshold}")
        
        print(f"\nList of contaminating genes saved to: {output_file}")

    except FileNotFoundError:
        print(f"Error: The file '{input_file}' was not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Filter contaminant genes from RNA-seq differential expression results.",
        formatter_class=argparse.RawTextHelpFormatter)
        
    parser.add_argument(
        '--input', 
        type=str, 
        required=True, 
        help="Path to the differential expression results file (CSV or TSV format).\nThis file should be the result of comparing the 'contaminated' group vs. the 'clean' group."
    )
    parser.add_argument(
        '--output', 
        type=str, 
        required=True, 
        help="Path for the output file, which will contain the list of contaminant gene IDs."
    )
    parser.add_argument(
        '--log2fc_threshold', 
        type=float, 
        default=5.0, 
        help="Log2 fold change threshold. Genes with a log2FC *above* this value will be flagged.\n(Default: 5.0)"
    )
    parser.add_argument(
        '--padj_threshold', 
        type=float, 
        default=0.01, 
        help="Adjusted p-value (padj) threshold. Genes with a padj *below* this value will be flagged.\n(Default: 0.01)"
    )

    args = parser.parse_args()

    filter_contaminants(args.input, args.output, args.log2fc_threshold, args.padj_threshold)

    print("\nTo run this script, use a command like this in your terminal:")
    print("python this_script_name.py --input your_de_results.csv --output contaminant_genes.txt --log2fc_threshold 5 --padj_threshold 0.01")