import pandas as pd

# Define the filtering thresholds
# A gene is a contaminant if its log2FC is less than this value
LOG2FC_THRESHOLD = -5.0  # Captures genes at least 32-fold higher in CAR_only samples
# And its average expression in the "clean" IL15 samples is below this value
COUNT_THRESHOLD = 20.0

# Load the DESeq2 results and normalized counts
results_df = pd.read_csv("initial_deseq_results.csv", index_col=0)
counts_df = pd.read_csv("normalized_counts.csv", index_col=0)

# Identify the sample columns for the CAR_IL15 group
il15_columns = [col for col in counts_df.columns if 'CAR_IL15' in col]

# Calculate the mean normalized count across the CAR_IL15 samples for each gene
counts_df['mean_il15_counts'] = counts_df[il15_columns].mean(axis=1)

# Apply the filtering logic to find contaminant genes
contaminant_mask = (results_df['log2FoldChange'] < LOG2FC_THRESHOLD) & \
                   (counts_df['mean_il15_counts'] < COUNT_THRESHOLD)

contaminant_genes = results_df[contaminant_mask]

print(f"Found {len(contaminant_genes)} potential contaminant genes.")

# Save the list of contaminant gene names to a file
contaminant_genes.index.to_series().to_csv("contaminant_genes.txt", index=False, header=False)

print("Contaminant gene list saved to contaminant_genes.txt")
