import pandas as pd
import numpy as np
import os

def filter_contaminants():
    """
    Filters a differential expression results file to remove contaminating
    cancer genes based on a log2FoldChange threshold derived from a known
    list of cancer genes.
    """
    # --- Step 1: Create dummy input files for demonstration ---
    # In your actual use, replace these sections with your real files.
    # Make sure your files are named 'de_results.tsv' and 'cancer_genes.txt'
    # and are in the same directory as the script.

    # Create a dummy de_results.tsv file
    de_data = {
        'gene_id': ['T_CELL_GENE_1', 'T_CELL_GENE_2', 'CANCER_GENE_1', 'CANCER_GENE_2', 
                    'BIOLOGICAL_DOWN_GENE', 'T_CELL_GENE_3', 'CANCER_GENE_3'],
        'log2FoldChange': [1.5, -0.5, -8.2, -9.5, -2.1, 2.5, -7.9],
        'pvalue': [0.001, 0.04, 0.00001, 0.000001, 0.03, 0.0001, 0.00002],
        'padj': [0.001, 0.04, 0.00001, 0.000001, 0.03, 0.0001, 0.00002]
    }
    pd.DataFrame(de_data).to_csv('de_results.tsv', sep='\t', index=False)

    # Create a dummy cancer_genes.txt file
    cancer_genes_list = ['CANCER_GENE_1', 'CANCER_GENE_2', 'CANCER_GENE_3']
    with open('cancer_genes.txt', 'w') as f:
        for gene in cancer_genes_list:
            f.write(gene + '\n')

    # --- Step 2: Load the data ---
    try:
        de_results_df = pd.read_csv('de_results.tsv', sep='\t')
        with open('cancer_genes.txt', 'r') as f:
            # Use a set for efficient lookup
            known_cancer_genes = {line.strip() for line in f if line.strip()}
    except FileNotFoundError as e:
        print(f"Error: Required file not found. Please ensure '{e.filename}' exists.")
        return

    # --- Step 3: Establish the LFC threshold from known contaminants ---
    # Find the rows in the DE results that correspond to our known cancer genes
    contaminant_df = de_results_df[de_results_df['gene_id'].isin(known_cancer_genes)]

    if contaminant_df.empty:
        print("No known cancer genes were found in the DE results file.")
        print("No filtering will be applied. Printing original results.")
        print(de_results_df.to_string())
        return

    # The threshold is the *maximum* LFC among the contaminants.
    # Any gene with an LFC this low or lower is suspect.
    lfc_threshold = contaminant_df['log2FoldChange'].max()

    # --- Step 4: Explain the filtering rule ---
    print("--- Contamination Filtering Report ---")
    print(f"Identified {len(contaminant_df)} known cancer genes in the DE results.")
    print(f"The maximum log2FoldChange among these contaminants is: {lfc_threshold:.4f}")
    print("\nThis value will be used as the filtering threshold.")
    print("The filtering rule is: Keep a gene if its 'log2FoldChange' is greater than the threshold.")
    print(f"Final Equation: Keep gene if log2FoldChange > {lfc_threshold:.4f}")
    print("--------------------------------------\n")

    # --- Step 5: Filter the data and print the output ---
    # The condition for keeping a gene is that its LFC must be > the threshold
    filtered_df = de_results_df[de_results_df['log2FoldChange'] > lfc_threshold].copy()
    
    print("Filtered Differential Expression Results:")
    # Print the DataFrame as a string, which is a clean way to show tabular data.
    print(filtered_df.to_string(index=False))

    # --- Step 6: Clean up dummy files ---
    os.remove('de_results.tsv')
    os.remove('cancer_genes.txt')

if __name__ == '__main__':
    filter_contaminants()