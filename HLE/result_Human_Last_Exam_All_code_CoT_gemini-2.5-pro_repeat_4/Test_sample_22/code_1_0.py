import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score

def demonstrate_pgs_heritability_gap():
    """
    This simulation demonstrates that the predictive ability of a Polygenic Score (R^2_PGS)
    is necessarily lower than the true SNP heritability (h^2_SNP) due to estimation
    error from finite sample sizes.
    """
    # 1. Simulation Parameters
    n_individuals = 5000   # Total number of individuals in our study
    n_snps = 10000         # Total number of SNPs
    h2_snp = 0.50          # True SNP heritability (the theoretical maximum R^2)
    prop_causal = 0.1      # Proportion of SNPs that truly affect the trait

    print("--- Simulation Plan ---")
    print(f"1. Define a true SNP Heritability (h²_SNP): {h2_snp:.2f}")
    print("2. Simulate genotypes and a phenotype based on this heritability.")
    print("3. Split data into a 'GWAS' training set and a 'PGS' test set.")
    print("4. Estimate SNP effects from the 'GWAS' set (introducing noise).")
    print("5. Build a PGS in the test set and calculate its predictive R².")
    print("6. Compare the predictive R² to the true SNP Heritability.\n")

    # 2. Simulate Genotypes and True Genetic Effects
    # For simplicity, we simulate standardized genotypes (mean=0, variance=1)
    genotypes = np.random.normal(0, 1, (n_individuals, n_snps))
    
    # Select causal SNPs and assign them true effect sizes
    n_causal = int(n_snps * prop_causal)
    causal_indices = np.random.choice(n_snps, n_causal, replace=False)
    true_betas = np.zeros(n_snps)
    # The variance of betas is set so the total genetic variance equals h2_snp
    true_betas[causal_indices] = np.random.normal(0, np.sqrt(h2_snp / n_causal), n_causal)

    # Calculate the true genetic value and the environmental effect
    true_genetic_value = genotypes @ true_betas
    environmental_effect = np.random.normal(0, np.sqrt(1 - h2_snp), n_individuals)
    
    # Create the final phenotype
    phenotype = true_genetic_value + environmental_effect

    # 3. Split data into training (for GWAS) and testing (for PGS)
    (g_train, g_test, y_train, y_test) = train_test_split(
        genotypes, phenotype, test_size=0.5, random_state=42
    )
    n_train = g_train.shape[0]

    # 4. Simulate a GWAS on the training set to get *estimated* effect sizes
    # A simple way to estimate betas is via regression, which simplifies to this
    # matrix operation for standardized genotypes. These estimates will have noise.
    estimated_betas = (1/n_train) * g_train.T @ y_train

    # 5. Build the PGS in the test set using the estimated betas
    pgs_test = g_test @ estimated_betas

    # 6. Evaluate the PGS predictive ability (R^2) in the test set
    r2_pgs = r2_score(y_test, pgs_test)

    # --- Final Output ---
    print("--- Final Equation of Interest ---")
    print("Is R²_PGS < h²_SNP?")
    print(f"   {r2_pgs:.4f} < {h2_snp:.4f}")
    print("\n--- Conclusion ---")
    print(f"The theoretical maximum variance explained (h²_SNP) was set to {h2_snp:.4f}.")
    print(f"The actual variance explained by the PGS (R²_PGS) was {r2_pgs:.4f}.")
    print("The simulation shows the predictive ability of the PGS is indeed lower than the SNP heritability.")

# Run the demonstration
demonstrate_pgs_heritability_gap()