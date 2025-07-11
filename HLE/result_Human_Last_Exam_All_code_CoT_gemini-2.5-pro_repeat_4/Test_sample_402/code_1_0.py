import pandas as pd
import io

# --- Step 1: Simulate the input data ---
# In a real scenario, you would load your data from a file like:
# df = pd.read_csv("differential_expression_results.csv")
#
# For this example, we'll create a sample dataset that mimics your results.
# This data includes:
# - Normal genes (e.g., GENE_A, GENE_B)
# - Upregulated T-cell genes in the IL15 group (e.g., T_CELL_GENE_1)
# - Potential contaminating cancer genes (e.g., CANCER_GENE_1, CANCER_GENE_2)
#   These have very large negative log2FoldChange values.
csv_data = """gene,log2FoldChange,padj
T_CELL_GENE_1,4.5,0.0001
GENE_A,0.5,0.6
CANCER_GENE_1,-9.2,0.00001
GENE_B,-1.1,0.04
T_CELL_GENE_2,3.1,0.002
GENE_C,0.1,0.9
CANCER_GENE_2,-7.8,0.00023
T_CELL_GENE_3,-2.5,0.03
"""

# Read the simulated data into a pandas DataFrame
df = pd.read_csv(io.StringIO(csv_data))

print("--- Initial Analysis ---")
print(f"Original number of genes: {len(df)}")
print("Original data:")
print(df.to_string())
print("\n" + "="*40 + "\n")

# --- Step 2: Define the filtering strategy and equation ---
# We will filter out genes that are highly specific to the control ("CAR only") group,
# which are presumed to be contaminants from the cancer cell line.
# The signature of these genes is a very large negative log2FoldChange and a significant p-value.
log2fc_threshold = -5.0
padj_threshold = 0.05

print("--- Filtering Process ---")
print("The filtering equation is: (log2FoldChange < log2fc_threshold) AND (padj < padj_threshold)")
print(f"Applying filter with values: (log2FoldChange < {log2fc_threshold}) AND (padj < {padj_threshold})\n")


# --- Step 3: Identify and report the contaminating genes ---
# Create a boolean mask to find genes that meet the contamination criteria.
is_contaminant = (df['log2FoldChange'] < log2fc_threshold) & (df['padj'] < padj_threshold)

contaminant_genes = df[is_contaminant]

if not contaminant_genes.empty:
    print("Identified potential contaminant genes to be removed:")
    print(contaminant_genes.to_string())
else:
    print("No genes met the contamination criteria.")

print("\n" + "="*40 + "\n")

# --- Step 4: Remove the contaminating genes and show the final result ---
# Use the inverse of the mask (~) to keep all genes that are NOT contaminants.
df_filtered = df[~is_contaminant]

print("--- Final Filtered Results ---")
print(f"Number of genes after filtering: {len(df_filtered)}")
print("Cleaned gene list (contaminants removed):")
print(df_filtered.to_string())
