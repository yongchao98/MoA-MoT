import numpy as np

def demonstrate_pgs_heritability():
    """
    This simulation demonstrates that the predictive ability of a Polygenic Score (R^2)
    is lower than the theoretical SNP heritability (h^2).
    """
    # --- 1. Simulation Parameters ---
    n_individuals = 5000
    n_snps = 1000
    # The proportion of variance in the trait explained by all SNPs (SNP heritability)
    snp_heritability = 0.50
    # Proportion of SNPs that are actually causal
    proportion_causal = 0.05

    # --- 2. Simulate Genotypes ---
    # Genotypes (0, 1, or 2 for allele counts) for each person at each SNP
    # We simulate based on a fictional allele frequency of 0.5
    genotypes = np.random.binomial(2, 0.5, size=(n_individuals, n_snps))

    # --- 3. Define True Genetic Architecture ---
    # Randomly select a subset of SNPs to be causal
    n_causal = int(n_snps * proportion_causal)
    causal_indices = np.random.choice(n_snps, n_causal, replace=False)

    # Assign true effect sizes (betas) to the causal SNPs
    true_betas = np.zeros(n_snps)
    true_betas[causal_indices] = np.random.normal(0, 1, size=n_causal)

    # --- 4. Calculate True Genetic Values and Phenotypes ---
    # Calculate the true genetic component of the phenotype for each individual
    true_genetic_value = genotypes @ true_betas

    # Calculate the environmental component needed to achieve the target heritability
    var_g = np.var(true_genetic_value)
    # Based on h^2 = Var(G) / Var(P), where Var(P) = Var(G) + Var(E)
    var_e = var_g * (1 - snp_heritability) / snp_heritability
    environmental_noise = np.random.normal(0, np.sqrt(var_e), size=n_individuals)

    # Create the final phenotype
    phenotype = true_genetic_value + environmental_noise

    # Verify our simulated heritability is close to the target
    empirical_h2 = var_g / np.var(phenotype)

    # --- 5. Simulate an Imperfect GWAS Result ---
    # In a real GWAS, beta estimates have noise. We simulate this by adding
    # noise to the true betas. This represents estimation error.
    noise_level = 0.8 # Higher noise simulates a smaller, less powerful GWAS
    estimation_noise = np.random.normal(0, noise_level * np.std(true_betas[causal_indices]), size=n_snps)
    estimated_betas = true_betas + estimation_noise

    # --- 6. Build the PGS and Calculate Its Predictive Power (R^2) ---
    # The PGS is built using the *estimated* (and therefore imperfect) betas
    pgs = genotypes @ estimated_betas

    # The predictive power (R^2) is the squared correlation between the PGS and the phenotype
    # This is equivalent to the proportion of variance in the phenotype explained by the PGS
    correlation = np.corrcoef(phenotype, pgs)[0, 1]
    pgs_r2 = correlation**2

    # --- 7. Print the Results ---
    print("--- Simulation Results ---")
    print(f"Target SNP Heritability (h^2): {snp_heritability:.4f}")
    print(f"Empirically Measured h^2 in simulation: {empirical_h2:.4f}\n")
    print("The SNP Heritability represents the THEORETICAL MAXIMUM variance that can be explained.")
    print("\nThe Polygenic Score (PGS) is a PRACTICAL model built from noisy estimates.")
    print(f"Predictive ability of the PGS (R^2): {pgs_r2:.4f}\n")
    print("Conclusion: The variance explained by the PGS is lower than the total SNP heritability.")
    print(f"{pgs_r2:.4f} < {snp_heritability:.4f}")


demonstrate_pgs_heritability()