import numpy as np

def demonstrate_pgs_vs_heritability():
    """
    A simulation to demonstrate that the variance explained by a Polygenic Score (PGS)
    is lower than the SNP heritability of a trait.
    """
    # --- 1. Simulation Parameters ---
    n_snps = 1000      # Number of SNPs in our simulated genome
    n_causal = 100     # Number of SNPs that actually affect the trait
    h2_snp = 0.50      # True SNP heritability (our theoretical maximum)
    n_gwas = 20000     # Sample size of the simulated discovery GWAS (to estimate effects)
    n_test = 2000      # Sample size of the testing cohort (to evaluate the PGS)
    
    # --- 2. Simulate a "True" Genetic Architecture ---
    # Create true effect sizes (betas) for all SNPs
    beta_true = np.zeros(n_snps)
    # Randomly select causal SNPs and assign them an effect size from a normal distribution
    causal_indices = np.random.choice(n_snps, n_causal, replace=False)
    beta_true[causal_indices] = np.random.normal(0, 1, n_causal)

    # --- 3. Simulate a Testing Cohort and their True Phenotypes ---
    # Simulate genotypes for the testing cohort (0, 1, or 2 copies of an allele)
    # We will standardize them to have mean=0 and variance=1 for simplicity
    maf = np.random.uniform(0.05, 0.5, n_snps)
    G_test = np.random.binomial(2, maf, size=(n_test, n_snps)).astype(float)
    G_test = (G_test - np.mean(G_test, axis=0)) / np.std(G_test, axis=0)

    # Calculate the "true" genetic component of the phenotype
    true_genetic_value = G_test @ beta_true
    
    # Scale true betas so the variance of the genetic component equals h2_snp
    scaling_factor = np.sqrt(h2_snp / np.var(true_genetic_value))
    beta_true_scaled = beta_true * scaling_factor
    true_genetic_value_scaled = G_test @ beta_true_scaled
    
    # Add environmental noise to create the final phenotype
    # The variance of environmental noise is set to 1 - h2_snp
    environmental_noise = np.random.normal(0, np.sqrt(1 - h2_snp), n_test)
    phenotype = true_genetic_value_scaled + environmental_noise

    # --- 4. Simulate GWAS Results (Imperfect Beta Estimates) ---
    # In a real GWAS, beta estimates have sampling error. Larger N = smaller error.
    # We simulate this by adding noise to the true betas.
    # The variance of the error is inversely proportional to GWAS sample size.
    error_variance = (1 - h2_snp) / n_gwas
    beta_error = np.random.normal(0, np.sqrt(error_variance), n_snps)
    beta_estimated = beta_true_scaled + beta_error
    
    # --- 5. Build the PGS and Measure its Predictive Ability ---
    # Calculate the PGS for the test cohort using the estimated betas
    pgs = G_test @ beta_estimated

    # Calculate the variance explained by the PGS (R-squared)
    # This is the squared correlation between the PGS and the true phenotype
    r_squared_pgs = np.corrcoef(pgs, phenotype)[0, 1] ** 2
    
    # --- 6. Print Results ---
    print(f"Theoretical Maximum (SNP Heritability, h²_snp): {h2_snp:.4f}")
    print(f"Achieved Prediction (PGS Variance Explained, R²): {r_squared_pgs:.4f}\n")
    print("As shown, the predictive ability of the PGS is lower than the SNP heritability.")

# Run the simulation and print the final answer
if __name__ == '__main__':
    demonstrate_pgs_vs_heritability()
