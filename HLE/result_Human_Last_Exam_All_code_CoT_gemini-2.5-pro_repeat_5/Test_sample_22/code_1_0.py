import numpy as np

def demonstrate_pgs_vs_heritability():
    """
    A simulation to demonstrate that the predictive ability of a Polygenic Score (PGS)
    is lower than the SNP heritability.
    """
    # --- 1. Simulation Parameters ---
    n_individuals = 5000
    n_snps = 1000
    n_causal_snps = 100
    true_heritability = 0.5  # This is the true SNP heritability (h²_SNP)

    # --- 2. Simulate a Population with a Defined Genetic Architecture ---
    print("Step 1: Simulating population genetics...")
    # Simulate genotypes (0, 1, or 2 copies of a minor allele)
    # Assume a minor allele frequency of 0.5 for simplicity
    genotypes = np.random.binomial(2, 0.5, size=(n_individuals, n_snps))
    # Standardize genotypes (mean 0, variance 1)
    genotypes = (genotypes - np.mean(genotypes, axis=0)) / np.std(genotypes, axis=0)

    # Define true causal SNPs and their effect sizes
    causal_indices = np.random.choice(n_snps, n_causal_snps, replace=False)
    true_effects = np.zeros(n_snps)
    # Effects drawn from a standard normal distribution
    true_effects[causal_indices] = np.random.randn(n_causal_snps)

    # Calculate the true genetic component (G) of the phenotype
    genetic_value = genotypes @ true_effects

    # Calculate the environmental component (E) to match the desired heritability
    genetic_variance = np.var(genetic_value)
    # Based on h² = Var(G) / (Var(G) + Var(E)), we derive Var(E)
    environmental_variance = genetic_variance * (1 / true_heritability - 1)
    environmental_noise = np.random.normal(0, np.sqrt(environmental_variance), n_individuals)

    # Create the final phenotype (P = G + E)
    phenotype = genetic_value + environmental_noise

    # Verify the simulated heritability (should be close to the target)
    simulated_heritability = genetic_variance / np.var(phenotype)
    print(f"   - Target SNP heritability (h²_SNP): {true_heritability:.4f}")
    print(f"   - Actual simulated h²_SNP: {simulated_heritability:.4f}\n")


    # --- 3. Simulate a GWAS and PGS Testing ---
    print("Step 2: Simulating GWAS and PGS testing...")
    # Split data into a discovery (GWAS) set and a testing set
    split_point = int(n_individuals * 0.8)
    gwas_geno, test_geno = genotypes[:split_point], genotypes[split_point:]
    gwas_pheno, test_pheno = phenotype[:split_point], phenotype[split_point:]

    # Run a "GWAS" on the discovery set to estimate SNP effects
    # (by regressing the phenotype on each SNP one at a time)
    estimated_effects = np.zeros(n_snps)
    for i in range(n_snps):
        # OLS estimate for single standardized regressor: cov(x,y)/var(x)
        # Since var(x) is 1, it's just cov(x,y)
        estimated_effects[i] = np.cov(gwas_geno[:, i], gwas_pheno)[0, 1]

    # Build the PGS in the testing set using the estimated effects
    pgs = test_geno @ estimated_effects

    # --- 4. Compare Predictive Ability vs. Heritability ---
    print("Step 3: Comparing results...")
    # Calculate the predictive ability (R²) of the PGS
    # R² is the squared correlation between predicted scores and true phenotypes
    corr_matrix = np.corrcoef(pgs, test_pheno)
    pgs_r_squared = corr_matrix[0, 1]**2

    print(f"\n--- Final Comparison ---")
    print(f"Theoretical Maximum (SNP Heritability): {simulated_heritability:.4f}")
    print(f"Achieved Prediction (PGS R²):         {pgs_r_squared:.4f}")
    print("\nAs shown, the predictive ability of the PGS is lower than the total SNP heritability.")

if __name__ == '__main__':
    demonstrate_pgs_vs_heritability()