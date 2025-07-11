import math

# --- Parameters ---
# SNP heritability (the theoretical maximum variance explained by SNPs)
h_snp_sq = 0.50

# Sample size of the GWAS used to derive the PGS weights
gwas_sample_size_N = 100000

# The effective number of independent SNPs contributing to the trait.
# This is often estimated to be around 60,000 for European ancestry.
effective_snps_M = 60000

# --- Calculation ---
# This formula predicts the expected R-squared of a polygenic score.
# R_sq_pgs = h_snp_sq / (1 + (M / (N * h_snp_sq)))

# Calculate the term in the denominator
denominator_term = gwas_sample_size_N * h_snp_sq
full_denominator = 1 + (effective_snps_M / denominator_term)

# Calculate the expected R-squared of the PGS
r_sq_pgs = h_snp_sq / full_denominator

# --- Output ---
print(f"This script calculates the expected predictive ability (R²) of a Polygenic Score (PGS).")
print("-" * 50)
print(f"Inputs:")
print(f"SNP Heritability (h²_SNP): {h_snp_sq}")
print(f"GWAS Sample Size (N): {gwas_sample_size_N}")
print(f"Effective Number of SNPs (M): {effective_snps_M}")
print("-" * 50)

print(f"The calculation is based on the formula:")
print(f"Expected R² = h²_SNP / (1 + M / (N * h²_SNP))")
print("\nSubstituting the numbers into the formula:")
# Here we print each number in the final equation
print(f"Expected R² = {h_snp_sq} / (1 + {effective_snps_M} / ({gwas_sample_size_N} * {h_snp_sq}))")
print("-" * 50)

print(f"Result:")
print(f"The calculated predictive ability of the PGS (R²) is: {r_sq_pgs:.4f}")
print(f"\nComparison:")
print(f"The PGS R² ({r_sq_pgs:.4f}) is lower than the SNP Heritability ({h_snp_sq:.4f}).")
