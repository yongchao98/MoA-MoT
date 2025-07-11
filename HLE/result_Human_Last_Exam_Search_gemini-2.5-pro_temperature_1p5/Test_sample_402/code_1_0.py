import pandas as pd
import argparse

# Set up command-line argument parsing
parser = argparse.ArgumentParser(description='Filter DGE results based on Log2 Fold Change.')
parser.add_argument('input_file', type=str, help='Path to the input DGE results file (CSV).')
parser.add_argument('output_file', type=str, help='Path for the filtered output file (CSV).')
parser.add_argument('--l2fc_column', type=str, default='log2FoldChange', help='Name of the L2FC column.')
parser.add_argument('--threshold', type=float, default=-5.0, help='The L2FC threshold.')

args = parser.parse_args()

# Read the DGE results file
try:
    df = pd.read_csv(args.input_file)
except FileNotFoundError:
    print(f"Error: Input file not found at {args.input_file}")
    exit()

# Apply the filter
filtered_df = df[df[args.l2fc_column] >= args.threshold]

# Save the filtered dataframe
filtered_df.to_csv(args.output_file, index=False)

print(f"Filtering complete. Original genes: {len(df)}. Filtered genes: {len(filtered_df)}.")
print(f"Results saved to {args.output_file}")