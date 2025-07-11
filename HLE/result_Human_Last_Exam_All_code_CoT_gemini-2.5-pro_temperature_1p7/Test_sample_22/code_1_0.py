import numpy as np
from sklearn.model_selection import train_test_split
from scipy.stats import pearsonr

def demonstrate_pgs_vs_heritability():
    """
    A simulation to demonstrate that the variance explained by a Polygenic Score (PGS)
    is lower than the SNP heritability of the trait.
    """
    # --- 1. Simulation Parameters ---
    n_individuals = 10000 # Number of individuals in our population
    n_snps = 2000         # Number of SNPs
    n_causal_snps = 200   # Number of SNPs that actually affect the trait
    h2_target = 0.60      # The true SNP heritability we are simulating

    print(f"Simulating with target SNP heritability: {h2_target:.2f}\n")

    # --- 2. Simulate Data ---
    # Simulate genotypes (0, 1, or 2 for allele count)
    genotypes = np.random.binomial(2, 0.3, size=(n_individuals, n_snps))

    # Simulate true SNP effect sizes (betas). Most are 0.
    true_betas = np.zeros(n_snps)
    causal_indices = np.random.choice(n_snps, n_causal_snps, replace=False)
    true_betas[causal_indices] = np.random.normal(0, 1, size=n_causal_snps)

    # Calculate the true genetic component of the trait
    genetic_value = genotypes @ true_betas

    # Calculate environmental component variance to meet the target heritability
    # h2 = Var(G) / (Var(G) + Var(E)) => Var(E) = Var(G) * (1 - h2) / h2
    var_g = np.var(genetic_value)
    var_e = var_g * (1 - h2_target) / h2_target if h2_target > 0 else np.inf

    # Create the environmental component and the final phenotype
    environmental_value = np.random.normal(0, np.sqrt(var_e), size=n_individuals)
    phenotype = genetic_value + environmental_value

    # --- 3. Calculate True SNP Heritability ---
    # This is the theoretical maximum R^2 possible
    snp_heritability = var_g / np.var(phenotype)

    # --- 4. Mimic GWAS and PGS Prediction ---
    # Split data into a "GWAS" training set and a "prediction" testing set
    X_train, X_test, y_train, y_test = train_test_split(
        genotypes, phenotype, test_size=0.5, random_state=42
    )

    # Estimate SNP effects from the "GWAS" training data (a simplified per-SNP regression)
    estimated_betas = np.zeros(n_snps)
    for i in range(n_snps):
        # Using Pearson correlation as a quick proxy for single-SNP effect size
        # This is a simplification but captures the spirit of a GWAS summary
        corr, _ = pearsonr(X_train[:, i], y_train)
        estimated_betas[i] = corr

    # Build the PGS in the testing set using the estimated effects
    pgs = X_test @ estimated_betas

    # --- 5. Calculate Predictive Ability of PGS ---
    # The predictive ability is the squared correlation between PGS and true phenotype
    # in the independent testing set.
    if np.var(pgs) == 0:
        pgs_r2 = 0
    else:
        pgs_r2 = pearsonr(pgs, y_test)[0] ** 2

    # --- 6. Final Comparison ---
    print("--- Results ---")
    print(f"Theoretical Maximum (SNP Heritability): {snp_heritability:.4f}")
    print(f"Achieved Prediction (PGS R-squared)  : {pgs_r2:.4f}")
    print("\nFinal comparison of the values shows:")
    print(f"PGS R-squared ({pgs_r2:.4f}) is lower than SNP Heritability ({snp_heritability:.4f})")
    print("\nTherefore, the statement is:")
    print(pgs_r2 < snp_heritability)


if __name__ == '__main__':
    demonstrate_pgs_vs_heritability()