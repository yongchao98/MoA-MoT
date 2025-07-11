import numpy as np

def simulate_pgs_performance():
    """
    A simulation to demonstrate that the predictive R^2 of a Polygenic Score (PGS)
    is lower than the SNP heritability (h_snp^2).
    """
    # --- 1. Simulation Parameters ---
    n_individuals_gwas = 5000  # Number of individuals in the discovery (GWAS) sample
    n_individuals_test = 1000  # Number of individuals in the target (testing) sample
    n_snps = 1000              # Number of SNPs
    h2_snp = 0.5               # True SNP heritability (50% of variance is from these SNPs)

    print(f"设定: 真实的SNP遗传力 (h²_snp) = {h2_snp:.2f}\n")

    # --- 2. Simulate Genotypes and True Effect Sizes ---
    # Simulate genotypes for both GWAS and test sets (0, 1, or 2 representing allele counts)
    # Assume a minor allele frequency of 0.25 for simplicity
    genotypes_gwas = np.random.binomial(2, 0.25, size=(n_individuals_gwas, n_snps))
    genotypes_test = np.random.binomial(2, 0.25, size=(n_individuals_test, n_snps))

    # Standardize genotypes (mean=0, variance=1)
    genotypes_gwas = (genotypes_gwas - np.mean(genotypes_gwas, axis=0)) / np.std(genotypes_gwas, axis=0)
    genotypes_test = (genotypes_test - np.mean(genotypes_test, axis=0)) / np.std(genotypes_test, axis=0)

    # Simulate true SNP effect sizes (betas)
    # The sum of squares of true betas for standardized genotypes equals the heritability
    true_betas = np.random.normal(0, np.sqrt(h2_snp / n_snps), n_snps)

    # --- 3. Simulate Phenotypes ---
    # Calculate the true genetic component of the phenotype
    true_genetic_value_gwas = genotypes_gwas @ true_betas
    true_genetic_value_test = genotypes_test @ true_betas

    # Simulate environmental noise (the remaining variance)
    environmental_variance = 1 - h2_snp
    noise_gwas = np.random.normal(0, np.sqrt(environmental_variance), n_individuals_gwas)
    noise_test = np.random.normal(0, np.sqrt(environmental_variance), n_individuals_test)

    # Create the final phenotype = genetic value + noise
    phenotype_gwas = true_genetic_value_gwas + noise_gwas
    phenotype_test = true_genetic_value_test + noise_test

    # --- 4. "Run" a GWAS to get Estimated Effect Sizes ---
    # In a real GWAS, we'd run a regression for each SNP.
    # A simplified way is to calculate the covariance, which is proportional to the regression coefficient
    # for standardized data. This mimics the estimation process with noise.
    # Estimated betas = true betas + estimation error.
    # The error decreases as GWAS sample size increases.
    estimated_betas = (genotypes_gwas.T @ phenotype_gwas) / n_individuals_gwas

    # --- 5. Build the PGS and Calculate its Predictive Ability (R^2) ---
    # Calculate the PGS for the test set using the *estimated* betas from the GWAS
    pgs_test = genotypes_test @ estimated_betas

    # Calculate the variance explained (R^2) by the PGS in the test set
    correlation = np.corrcoef(pgs_test, phenotype_test)[0, 1]
    pgs_r2 = correlation**2

    print("--- 模拟结果 ---")
    print(f"理论上的最大可解释方差 (SNP遗传力): {h2_snp:.4f}")
    print(f"基于模拟GWAS数据构建的多基因评分 (PGS) 的实际预测方差 (R²): {pgs_r2:.4f}\n")
    print("结论: 正如预期的那样，由于GWAS样本量有限导致效应量估计不完美，PGS的预测能力 (R²) 低于理论上的SNP遗传力。")

if __name__ == "__main__":
    simulate_pgs_performance()