import pandas as pd
import argparse
import sys

def filter_contaminating_genes(input_file, output_file, lfc_threshold, padj_threshold):
    """
    Filters contaminating genes from a differential expression results file.

    Contaminants are defined as genes with a log2FoldChange < lfc_threshold
    and padj < padj_threshold.
    """
    try:
        # Read the differential expression data
        # The script assumes a CSV/TSV file with columns 'log2FoldChange' and 'padj'
        # It will try to auto-detect the separator (comma or tab)
        try:
            df = pd.read_csv(input_file)
            if df.shape[1] == 1:
                df = pd.read_csv(input_file, sep='\t')
        except Exception as e:
            print(f"Error reading file: {e}", file=sys.stderr)
            print("Please ensure the input file is a valid CSV or TSV.", file=sys.stderr)
            sys.exit(1)

        # Ensure required columns exist
        required_cols = ['log2FoldChange', 'padj']
        if not all(col in df.columns for col in required_cols):
            print(f"Error: Input file must contain the columns: {required_cols}", file=sys.stderr)
            sys.exit(1)
            
        # --- Filtering Logic ---
        # The comparison is assumed to be (CAR T + IL15) vs (CAR T only).
        # Contaminating genes from cancer cells in the "CAR T only" samples will have
        # a large negative log2FoldChange.
        
        print("--- Contamination Filtering Strategy ---")
        print(f"Analysis assumes comparison: (CAR T + IL15) vs (CAR T only)")
        print("Contaminating genes are identified by high expression in the 'CAR T only' group.")
        print("\nThe filtering equation to identify a contaminating gene is:")
        print(f"log2FoldChange < {lfc_threshold} AND padj < {padj_threshold}")
        
        # Create a boolean mask to identify the contaminating genes
        contaminant_mask = (df['log2FoldChange'] < lfc_threshold) & (df['padj'] < padj_threshold)
        
        contaminant_genes = df[contaminant_mask]
        
        if contaminant_genes.empty:
            print("\nNo contaminating genes found with the given thresholds.")
        else:
            print(f"\nFound {len(contaminant_genes)} potential contaminating genes to remove:")
            # Assuming a column with gene names/IDs exists, print them.
            # Trying common names for the gene identifier column.
            gene_col_options = ['gene_name', 'gene', 'Gene', 'id', 'ID']
            gene_col = next((col for col in gene_col_options if col in df.columns), df.columns[0])
            for gene in contaminant_genes[gene_col]:
                print(f"- {gene}")

        # Filter the dataframe to exclude the contaminating genes
        filtered_df = df[~contaminant_mask]
        
        # Save the cleaned dataframe to the output file
        filtered_df.to_csv(output_file, index=False)
        
        print("\n--- Summary ---")
        print(f"Original number of genes: {len(df)}")
        print(f"Number of genes removed: {len(contaminant_genes)}")
        print(f"Final number of genes: {len(filtered_df)}")
        print(f"\nSuccessfully saved the filtered data to: {output_file}")

    except FileNotFoundError:
        print(f"Error: The file '{input_file}' was not found.", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Filter contaminating genes from RNA-seq differential expression results. "
                    "This script assumes the comparison was performed as (Treated vs Control), "
                    "where the Control group is the one with the contamination.",
        formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument(
        '-i', '--input', 
        required=True, 
        help="Path to the input differential expression file (CSV or TSV)."
    )
    parser.add_argument(
        '-o', '--output', 
        required=True, 
        help="Path for the output filtered file (CSV)."
    )
    parser.add_argument(
        '--lfc', 
        type=float, 
        default=-5.0,
        help="Log2 fold change threshold. Genes with a log2FC *below* this value will be considered. \nDefault: -5.0"
    )
    parser.add_argument(
        '--padj', 
        type=float, 
        default=0.05,
        help="Adjusted p-value (FDR) threshold. Genes with a padj *below* this value will be considered. \nDefault: 0.05"
    )
    
    args = parser.parse_args()
    
    filter_contaminating_genes(args.input, args.output, args.lfc, args.padj)
