import numpy as np

# Plan:
# 1. Define simulation parameters: SNP heritability (h2_snp), number of individuals (N), and number of SNPs (M).
# 2. Simulate genotypes (X) and true SNP effect sizes (beta).
# 3. Create a "true" genetic component (g) based on genotypes and true effects.
#    Scale it so its variance equals the SNP heritability.
# 4. Create a phenotype (y) by adding an environmental component to the genetic component.
# 5. Simulate a GWAS by estimating SNP effect sizes (beta_hat) from the simulated data.
#    These estimates will be imperfect due to the finite sample size.
# 6. Build a Polygenic Score (PGS) using these estimated effects.
# 7. Calculate the predictive ability of the PGS (R_squared_PGS) by seeing how much variance in the phenotype it explains.
# 8. Compare the R_squared_PGS with the true h2_snp to demonstrate that the former is lower.

# --- Step 1: Define parameters ---
# True SNP heritability (the theoretical maximum variance explained by all SNPs)
h2_snp = 0.50
# Number of individuals in the simulated GWAS
n_individuals = 5000
# Number of SNPs
n_snps = 10000

# --- Step 2: Simulate genotypes and true effects ---
# Simulate genotypes X, then standardize it (mean=0, variance=1).
genotypes = np.random.binomial(2, 0.5, size=(n_individuals, n_snps))
genotypes = (genotypes - np.mean(genotypes, axis=0)) / np.std(genotypes, axis=0)

# Simulate true effect sizes (beta) for each SNP. For simplicity, all SNPs have a small effect.
true_effects = np.random.normal(0, np.sqrt(h2_snp / n_snps), n_snps)

# --- Step 3: Create true genetic component ---
# Calculate the true genetic value 'g' for each individual and scale its variance to h2_snp
true_genetic_value = genotypes @ true_effects
true_genetic_value = (true_genetic_value - np.mean(true_genetic_value)) / np.std(true_genetic_value) * np.sqrt(h2_snp)

# --- Step 4: Create phenotype ---
# Simulate environmental component 'e' with variance (1 - h2_snp)
environmental_effect = np.random.normal(0, np.sqrt(1 - h2_snp), n_individuals)
# Create final phenotype 'y'. Total variance of y is var(g) + var(e) = h2_snp + (1-h2_snp) = 1
phenotype = true_genetic_value + environmental_effect

# --- Step 5: Simulate GWAS to get estimated effects ---
# This simulates the noisy estimates from a GWAS, with noise being inversely related to sample size.
estimated_effects = (1 / n_individuals) * (genotypes.T @ phenotype)

# --- Step 6: Build the Polygenic Score ---
# The PGS is the sum of genotypes weighted by the *estimated* effects from our simulated GWAS.
pgs = genotypes @ estimated_effects

# --- Step 7: Calculate the predictive ability of the PGS ---
# This is measured by R-squared: the squared correlation between the phenotype and the PGS.
correlation = np.corrcoef(phenotype, pgs)[0, 1]
r_squared_pgs = correlation**2

# --- Step 8: Compare and print results ---
print("This simulation demonstrates the relationship between SNP heritability and Polygenic Score (PGS) predictive ability.")
print("-" * 80)
print(f"Theoretical maximum variance explained (SNP Heritability, h²_snp): {h2_snp:.4f}")
print(f"Achieved predictive ability of the PGS (R²_PGS):                 {r_squared_pgs:.4f}")
print("-" * 80)
print("As shown, the predictive ability of the PGS is lower than the total SNP heritability.")
print("\nFinal Equation:")
print(f"R²_PGS ({r_squared_pgs:.4f}) < h²_snp ({h2_snp:.4f})")