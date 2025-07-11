import numpy as np

def demonstrate_pgs_heritability():
    """
    A simulation to demonstrate that the predictive R^2 of a Polygenic Score (PGS)
    is lower than the SNP heritability.
    """
    # --- 1. Simulation Setup ---
    # We define parameters for our simulated population and trait.
    n_individuals = 2000  # Number of people in our study
    n_snps = 5000       # Number of SNPs we measure
    h2_snp = 0.60       # The "true" SNP heritability we are aiming for
    
    # --- 2. Simulate Genotypes and True Genetic Effects ---
    # Create random genotypes (0, 1, or 2) for each person at each SNP
    genotypes = np.random.randint(0, 3, size=(n_individuals, n_snps))
    # Standardize genotypes for modeling (mean=0, variance=1)
    genotypes = (genotypes - np.mean(genotypes, axis=0)) / np.std(genotypes, axis=0)

    # Define a small fraction of SNPs as having a true causal effect
    n_causal_snps = int(n_snps * 0.01) # 1% of SNPs are causal
    true_betas = np.zeros(n_snps)
    causal_indices = np.random.choice(n_snps, n_causal_snps, replace=False)
    # The true effect sizes are drawn from a normal distribution
    true_betas[causal_indices] = np.random.normal(0, 1, n_causal_snps)

    # --- 3. Simulate Phenotype based on True Genetics + Environment ---
    # The true genetic value (G) is the sum of genotypes weighted by true effects
    G = genotypes @ true_betas
    var_G = np.var(G)

    # The environmental noise (E) makes up the rest of the variance
    # h2 = Var(G) / (Var(G) + Var(E))  =>  Var(E) = Var(G) * (1/h2 - 1)
    var_E = var_G * (1 / h2_snp - 1)
    E = np.random.normal(0, np.sqrt(var_E), n_individuals)
    
    # The final phenotype (P) is the sum of genetic and environmental components
    phenotype = G + E
    var_P = np.var(phenotype)

    # --- 4. Calculate the "True" SNP Heritability ---
    # This is the theoretical maximum variance we could ever explain with these SNPs.
    calculated_h2_snp = var_G / var_P

    # --- 5. Simulate a "GWAS" and build a PGS ---
    # In a real GWAS, effect sizes are estimated with error. We simulate this
    # by adding noise to the "true" betas.
    # The amount of noise is inversely related to GWAS power/sample size.
    error_variance = 0.8
    gwas_error = np.random.normal(0, error_variance, n_snps)
    estimated_betas = true_betas + gwas_error
    
    # Build the polygenic score (PGS) using these imperfect, estimated betas
    pgs = genotypes @ estimated_betas

    # --- 6. Calculate the Predictive Ability (R^2) of the PGS ---
    # The R^2 is the squared correlation between the actual phenotype and our PGS.
    correlation = np.corrcoef(phenotype, pgs)[0, 1]
    pgs_r2 = correlation**2

    # --- 7. Print the Comparison ---
    print("--- Simulation Results ---")
    print("\nStep 1: Calculating the true SNP Heritability (h²_SNP)")
    print(f"This is the theoretical maximum predictive power.")
    print(f"h²_SNP = Variance(Genetics) / Variance(Phenotype)")
    print(f"h²_SNP = {var_G:.4f} / {var_P:.4f} = {calculated_h2_snp:.4f}")
    
    print("\nStep 2: Calculating the predictive ability of the Polygenic Score (R²)")
    print("This is the actual predictive power achieved using estimated effects.")
    print(f"R² = (Correlation(Phenotype, PGS))²")
    print(f"R² = ({correlation:.4f})² = {pgs_r2:.4f}")

    print("\n--- Conclusion ---")
    print(f"The variance explained by the PGS (R² = {pgs_r2:.4f}) is lower than the total SNP heritability (h²_SNP = {calculated_h2_snp:.4f}).")
    print("This demonstrates that a PGS, built on imperfect estimates, captures a fraction of the total heritable genetic signal.")

if __name__ == '__main__':
    demonstrate_pgs_heritability()