import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split
import scipy.stats

def simulate_pgs_heritability_comparison():
    """
    This function simulates genetic data to demonstrate that the predictive
    ability of a Polygenic Score (PGS) is lower than the SNP heritability.
    """
    # 1. Set simulation parameters
    n_individuals = 5000
    n_snps = 1000
    n_causal_snps = 100
    target_h2 = 0.5  # We are aiming for a true heritability of 50%

    # 2. Simulate genotypes (coded as 0, 1, 2)
    # Using random minor allele frequencies (MAFs) for each SNP
    mafs = np.random.uniform(0.05, 0.5, n_snps)
    genotypes = np.random.binomial(2, mafs, size=(n_individuals, n_snps)).astype(float)

    # 3. Simulate true genetic effects and the "true genetic value"
    true_betas = np.zeros(n_snps)
    causal_indices = np.random.choice(n_snps, n_causal_snps, replace=False)
    # Effect sizes drawn from a standard normal distribution for causal SNPs
    true_betas[causal_indices] = np.random.normal(0, 1, n_causal_snps)
    
    # Calculate the true genetic value and standardize it to have variance of 1
    genetic_value = genotypes @ true_betas
    genetic_value = (genetic_value - np.mean(genetic_value)) / np.std(genetic_value)

    # 4. Simulate environmental noise to achieve the target heritability
    # Var(Phenotype) = Var(Genetic) + Var(Environment)
    # h2 = Var(Genetic) / (Var(Genetic) + Var(Environment))
    # Since Var(Genetic) is 1, h2 = 1 / (1 + Var(Environment))
    # Var(Environment) = (1 / h2) - 1
    var_env = (1 / target_h2) - 1
    std_env = np.sqrt(var_env)
    environmental_noise = np.random.normal(0, std_env, n_individuals)
    
    # 5. Create final phenotype
    phenotype = genetic_value + environmental_noise

    # 6. Calculate the true SNP heritability in our simulated population
    # This is the theoretical maximum variance that can be explained by our SNPs.
    snp_heritability = np.var(genetic_value) / np.var(phenotype)

    # 7. Split data into a discovery (GWAS) set and a target (prediction) set
    gwas_geno, target_geno, gwas_pheno, target_pheno = train_test_split(
        genotypes, phenotype, test_size=0.5, random_state=42)

    # 8. Simulate a GWAS by estimating effect sizes in the discovery set
    # Note: This is a simplified per-SNP regression
    estimated_betas = np.zeros(n_snps)
    for i in range(n_snps):
        snp_column = gwas_geno[:, i].reshape(-1, 1)
        model = LinearRegression().fit(snp_column, gwas_pheno)
        estimated_betas[i] = model.coef_[0]

    # 9. Build the PGS in the target set using the estimated effects
    pgs = target_geno @ estimated_betas

    # 10. Calculate the predictive ability (R^2) of the PGS
    # R^2 is the squared correlation between the PGS and the phenotype in the target set.
    r_val, _ = scipy.stats.pearsonr(pgs, target_pheno)
    pgs_r_squared = r_val**2
    
    # 11. Print the results for comparison
    print("--- Simulation Results ---")
    print(f"Theoretical maximum (SNP Heritability): h²_SNP = {snp_heritability:.4f}")
    print(f"Achieved prediction (PGS R-squared):     R²_PGS = {pgs_r_squared:.4f}")
    print("\nConclusion from simulation:")
    print("The predictive ability of the PGS (R²_PGS) is lower than the SNP heritability (h²_SNP).")
    print("This is because the PGS is built using *estimated* effect sizes from a finite GWAS sample, which contain statistical noise.")
    print("The SNP heritability represents the theoretical ceiling based on *true* effect sizes.")


if __name__ == '__main__':
    simulate_pgs_heritability_comparison()
