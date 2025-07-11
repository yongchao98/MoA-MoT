import numpy as np
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score

def simulate_pgs_vs_heritability():
    """
    A simulation to demonstrate that the variance explained by a Polygenic Score (PGS R^2)
    is lower than the true SNP heritability (h^2_SNP).
    """
    # --- 1. Simulation Parameters ---
    # The true proportion of variance explained by all SNPs (SNP Heritability)
    h2_snp = 0.50

    # Number of SNPs
    n_snps = 1000

    # Sample sizes for our simulated "studies"
    n_discovery = 5000  # Cohort to estimate SNP effects (our "GWAS")
    n_target = 2000     # Cohort to test the PGS

    print(f"Starting simulation with the following parameters:")
    print(f"True SNP Heritability (h²_SNP): {h2_snp}")
    print(f"Number of SNPs: {n_snps}")
    print(f"Discovery sample size: {n_discovery}")
    print(f"Target sample size: {n_target}")
    print("-" * 30)

    # --- 2. Generate "True" SNP Effects ---
    # In nature, we don't know these. But in a simulation, we can define them.
    # We draw them from a normal distribution. The scaling ensures the total genetic
    # variance adds up to h2_snp.
    np.random.seed(42)
    true_effects = np.random.normal(0, np.sqrt(h2_snp / n_snps), n_snps)

    # --- 3. Create Discovery Cohort and "Run GWAS" ---
    # Generate random genotypes (0, 1, or 2 copies of the minor allele)
    maf = np.random.uniform(0.05, 0.5, n_snps)
    G_discovery = np.hstack([np.random.binomial(2, p, (n_discovery, 1)) for p in maf])
    G_discovery_std = (G_discovery - G_discovery.mean(axis=0)) / G_discovery.std(axis=0)

    # Generate phenotypes based on true genetics and some random environmental noise
    genetic_component_discovery = G_discovery_std @ true_effects
    environmental_noise_discovery = np.random.normal(0, np.sqrt(1 - h2_snp), n_discovery)
    y_discovery = genetic_component_discovery + environmental_noise_discovery

    # "Run GWAS": Use a linear model to *estimate* the SNP effects from discovery data.
    # We use Ridge regression, which is common for this type of problem.
    gwas_model = Ridge(alpha=1.0)
    gwas_model.fit(G_discovery_std, y_discovery)
    estimated_effects = gwas_model.coef_
    print("GWAS simulation complete. SNP effects have been estimated.")

    # --- 4. Create Target Cohort and Calculate PGS R^2 ---
    # Generate a new, independent set of people
    G_target = np.hstack([np.random.binomial(2, p, (n_target, 1)) for p in maf])
    G_target_std = (G_target - G_target.mean(axis=0)) / G_target.std(axis=0)

    # Generate their true phenotypes
    genetic_component_target = G_target_std @ true_effects
    environmental_noise_target = np.random.normal(0, np.sqrt(1 - h2_snp), n_target)
    y_target = genetic_component_target + environmental_noise_target

    # Calculate the Polygenic Score for the target cohort using the *estimated* effects
    pgs_target = G_target_std @ estimated_effects
    print("Polygenic Score (PGS) calculated for the target cohort.")

    # Measure the predictive ability of the PGS (its R^2)
    pgs_r2 = r2_score(y_target, pgs_target)
    print("-" * 30)

    # --- 5. Compare and Conclude ---
    print("Final Comparison:")
    # The "equation" showing the relationship between the two values
    print(f"Predictive Ability (PGS R²)  <  Theoretical Maximum (SNP Heritability)")
    print(f"{pgs_r2:.4f}  <  {h2_snp:.4f}")

    if pgs_r2 < h2_snp:
        print("\nConclusion: As expected, the variance explained by the PGS is lower than the true SNP heritability.")
    else:
        print("\nConclusion: The simulation produced an unexpected result.")

if __name__ == '__main__':
    simulate_pgs_vs_heritability()