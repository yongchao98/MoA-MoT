import pandas as pd
import io

# --- Step 1: Prepare Sample Data ---
# In a real-world scenario, you would load your data from a file like this:
# df = pd.read_csv("path/to/your/de_results.csv")
#
# For this example, we'll create a sample dataset in memory that simulates
# your differential expression results. It includes some T-cell-related
# genes (CD247, GZMB, PRF1) and two likely cancer-cell contaminants.
csv_data = """gene_name,baseMean,log2FoldChange,pvalue,padj
GENE_A,150.5,2.5,0.001,0.01
CD247,450.2,1.5,0.02,0.04
CANCER_GENE_1,1000.1,-8.2,0.0001,0.001
GZMB,75.9,3.1,0.0005,0.008
CANCER_GENE_2,500.7,-10.5,0.00001,0.0002
PRF1,810.3,2.9,0.0004,0.007
"""

# Use io.StringIO to simulate reading from a file
data_file = io.StringIO(csv_data)
df = pd.read_csv(data_file)

# --- Step 2: Define the Filtering Strategy ---
# We will filter out genes with a log2 Fold Change below a certain negative threshold.
# This indicates high expression in the control (denominator) and low expression
# in the treatment (numerator), which matches the contamination pattern.
log2fc_threshold = -5.0

# --- Step 3: Apply the Filter ---
# The filtering 'equation' is based on the log2FoldChange column.
# First, identify the genes that meet this contamination criteria.
contaminant_genes = df[df['log2FoldChange'] < log2fc_threshold]

# Next, create a new dataframe that excludes these contaminating genes.
# The '~' operator inverts the boolean selection, keeping everything that
# does NOT meet the contamination criteria.
cleaned_df = df[~(df['log2FoldChange'] < log2fc_threshold)]


# --- Step 4: Output the Results ---
print("--- Filtering Strategy ---")
print("Contaminating genes are identified as those highly expressed in the control group (CAR T)")
print("and lowly expressed in the treatment group (CAR T + IL15). This corresponds to a large")
print("negative log2 Fold Change.")
print("\n--- Filtering Equation ---")
# As requested, printing the filtering logic and the number used.
print(f"Genes to be removed are selected where: log2FoldChange < {log2fc_threshold}")

print("\n--- Identified Contaminant Genes (Removed) ---")
if not contaminant_genes.empty:
    # We use .to_string() to print the column without the index.
    print(contaminant_genes['gene_name'].to_string(index=False))
else:
    print("No contaminant genes found with the current threshold.")

print("\n--- Cleaned Differential Expression Data ---")
# Print the final, filtered dataframe to the console.
# You can redirect this output to a new file using '>' in the command line.
# e.g., python your_script.py > cleaned_de_results.csv
print(cleaned_df.to_string(index=False))