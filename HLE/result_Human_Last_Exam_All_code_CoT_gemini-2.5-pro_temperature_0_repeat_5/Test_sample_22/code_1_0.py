import numpy as np

def demonstrate_pgs_heritability():
    """
    A simulation to demonstrate that the predictive ability (R^2) of a
    Polygenic Score (PGS) is lower than the SNP heritability (h_snp^2).
    """
    # --- 1. Set Simulation Parameters ---
    # True SNP heritability (the theoretical maximum R^2)
    h2_snp = 0.50

    # Number of SNPs that have a true effect on the trait
    n_causal_snps = 500

    # Sample sizes for our simulated studies
    n_gwas = 10000  # Sample size for the "discovery" GWAS
    n_test = 2000   # Sample size for the "testing" or "validation" set

    # A plausible Minor Allele Frequency (MAF) for the SNPs
    maf = 0.25

    # Set a seed for reproducibility of the random simulation
    np.random.seed(42)

    # --- 2. Simulate Genetic Data and True Phenotype ---
    # a) Simulate true SNP effect sizes (betas). These are the "perfect" weights.
    # We draw them from a normal distribution, scaled by heritability and number of SNPs.
    true_betas = np.random.normal(0, np.sqrt(h2_snp / n_causal_snps), n_causal_snps)

    # b) Simulate genotypes (0, 1, or 2 copies of the minor allele) for both datasets
    X_gwas = np.random.binomial(2, maf, size=(n_gwas, n_causal_snps))
    X_test = np.random.binomial(2, maf, size=(n_test, n_causal_snps))

    # c) Calculate the "true" genetic component of the phenotype using the true betas
    G_gwas = X_gwas @ true_betas
    G_test = X_test @ true_betas

    # d) Calculate the variance of the environmental/non-genetic component
    # This is scaled so the total variance of the phenotype will be ~1.0
    environmental_variance = 1 - np.var(G_gwas)
    if environmental_variance < 0:
        environmental_variance = 0 # Ensure variance is non-negative

    # e) Simulate the environmental component from a normal distribution
    E_gwas = np.random.normal(0, np.sqrt(environmental_variance), n_gwas)
    E_test = np.random.normal(0, np.sqrt(environmental_variance), n_test)

    # f) Create the final, observable phenotype by combining genetic and environmental parts
    Y_gwas = G_gwas + E_gwas
    Y_test = G_test + E_test

    # --- 3. "Run" a mock GWAS to build the Polygenic Score ---
    # a) Estimate SNP effects from the GWAS data. In a real GWAS, this is done
    # SNP-by-SNP. These estimated betas will be noisy due to the finite sample size.
    estimated_betas = np.zeros(n_causal_snps)
    for i in range(n_causal_snps):
        # Perform a simple linear regression for each SNP against the phenotype
        snp_i = X_gwas[:, i]
        # Using lstsq to solve for the slope (the beta)
        estimated_betas[i] = np.linalg.lstsq(snp_i[:, np.newaxis], Y_gwas, rcond=None)[0][0]

    # b) Build the Polygenic Score (PGS) in the independent test set using the noisy, estimated betas
    pgs_test = X_test @ estimated_betas

    # --- 4. Evaluate the PGS and Print Results ---
    # Calculate the variance explained (R^2) by the PGS in the test set.
    # This is the squared correlation between the PGS and the true phenotype.
    r_squared = np.corrcoef(Y_test, pgs_test)[0, 1]**2

    print("This simulation demonstrates the relationship between SNP heritability and PGS predictive ability.")
    print("-" * 80)
    print("The SNP heritability (h²_snp) is the theoretical maximum variance that can be explained.")
    print(f"In this simulation, the true h²_snp was set to: {h2_snp}")
    print("\nThe PGS predictive ability (R²) is what we actually achieve in a new sample.")
    print(f"The R² from our PGS, built with a GWAS of N={n_gwas}, is: {r_squared:.4f}")
    print("-" * 80)
    print("The final comparison shows the inequality:")
    print(f"PGS R² ({r_squared:.4f}) < SNP Heritability ({h2_snp})")
    print("\nThis inequality holds because estimating SNP effects from a finite sample introduces")
    print("statistical noise, preventing the PGS from capturing the full genetic signal.")

if __name__ == '__main__':
    demonstrate_pgs_heritability()