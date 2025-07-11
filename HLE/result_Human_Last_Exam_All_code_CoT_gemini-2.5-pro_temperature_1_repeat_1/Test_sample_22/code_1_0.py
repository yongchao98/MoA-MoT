import numpy as np
from scipy.stats import pearsonr

def demonstrate_pgs_heritability():
    """
    A simulation to demonstrate that the predictive ability of a Polygenic Score (PGS)
    is lower than the SNP heritability.
    """
    # 1. Simulation parameters
    n_individuals_gwas = 20000  # Sample size for our simulated "GWAS"
    n_individuals_test = 2000   # Sample size for testing the PGS
    n_snps = 5000               # Total number of SNPs in the genome
    n_causal_snps = 500         # Number of SNPs that truly affect the phenotype

    # 2. Generate genetic data
    # Create random allele frequencies for all SNPs
    mafs = np.random.uniform(0.05, 0.5, n_snps)
    # Generate genotypes (0, 1, or 2 copies of an allele) for the test set
    genotypes_test = np.random.binomial(2, mafs, size=(n_individuals_test, n_snps))
    # Standardize genotypes to have a mean of 0 and standard deviation of 1
    genotypes_test_std = (genotypes_test - np.mean(genotypes_test, axis=0)) / np.std(genotypes_test, axis=0)

    # 3. Define the "true" genetic model and simulate phenotype
    # Randomly choose which SNPs are causal
    causal_indices = np.random.choice(n_snps, n_causal_snps, replace=False)
    # Define true effect sizes (betas) for all SNPs (zero for non-causal ones)
    true_betas = np.zeros(n_snps)
    # Causal effects are drawn from a normal distribution
    true_betas[causal_indices] = np.random.normal(0, np.sqrt(1/n_causal_snps), n_causal_snps)

    # Calculate the true genetic component of the phenotype
    true_genetic_value = genotypes_test_std @ true_betas
    genetic_variance = np.var(true_genetic_value)

    # We set environmental variance to be equal to genetic variance,
    # targeting a true heritability of 0.5
    environmental_variance = genetic_variance
    environmental_noise = np.random.normal(0, np.sqrt(environmental_variance), n_individuals_test)

    # The final phenotype is the sum of genetic and environmental effects
    phenotype = true_genetic_value + environmental_noise
    total_variance = np.var(phenotype)

    # SNP heritability is the proportion of total variance explained by genetics
    snp_heritability = genetic_variance / total_variance

    # 4. Simulate a GWAS to get imperfect effect sizes for the PGS
    # The error in beta estimates is inversely related to the GWAS sample size
    error_variance = (total_variance - genetic_variance) / n_individuals_gwas
    estimation_error = np.random.normal(0, np.sqrt(error_variance), n_snps)
    estimated_betas = true_betas + estimation_error

    # 5. Build the PGS and calculate its predictive ability (R-squared)
    # The PGS is the sum of genotypes weighted by the *estimated* effects
    pgs = genotypes_test_std @ estimated_betas
    # R-squared is the squared correlation between the true phenotype and the PGS
    r_squared_pgs = pearsonr(phenotype, pgs)[0]**2

    # 6. Print the results
    print("This simulation compares the theoretical SNP heritability with the practical predictive ability of a Polygenic Score (PGS).")
    print("-" * 70)
    print(f"Theoretical maximum (SNP heritability, h²_SNP): {snp_heritability:.4f}")
    print(f"Achieved prediction (PGS R-squared, R²_PGS):     {r_squared_pgs:.4f}")
    print("-" * 70)
    print("As shown, the predictive R² of the PGS is lower than the total SNP heritability.")
    print("This is because the PGS is built using imperfect effect sizes from a finite-sized GWAS.\n")
    print("Final equation:")
    print(f"{r_squared_pgs:.4f} (R²_PGS) < {snp_heritability:.4f} (h²_SNP)")

if __name__ == '__main__':
    demonstrate_pgs_heritability()