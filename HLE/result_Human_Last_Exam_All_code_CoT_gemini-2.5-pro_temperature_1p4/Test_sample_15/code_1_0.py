import numpy as np

def run_simulation():
    """
    Simulates the described bioinformatics scenario to test for bias in population genetics estimators.
    """
    # --- 1. Simulation Parameters ---
    n_samples = 50  # Number of diploid individuals
    n_sites = 10000 # Length of the sequence
    true_polymorphic_sites = 100 # Number of true SNVs in the population
    p_filter = 0.3  # Probability of filtering a true alternative allele due to "low quality"

    n_haplotypes = 2 * n_samples

    print(f"--- Simulation Setup ---")
    print(f"Simulating {n_samples} samples ({n_haplotypes} haplotypes) over {n_sites} sites.")
    print(f"Introducing {true_polymorphic_sites} true polymorphic sites (SNVs).")
    print(f"Each alternative allele has a {p_filter*100:.1f}% chance of being filtered and imputed with the reference.\n")

    # --- 2. Generate the "True" Haplotype Data ---
    # Start with all reference alleles (represented by 0)
    true_haplotypes = np.zeros((n_haplotypes, n_sites), dtype=int)

    # Introduce polymorphic sites (alternative allele represented by 1)
    # Ensure all have a minor allele frequency > 0
    segregating_site_indices = np.random.choice(range(n_sites), true_polymorphic_sites, replace=False)
    for site in segregating_site_indices:
        # Assign alt alleles to a random number of haplotypes (between 1 and n_haplotypes/2)
        num_alt_alleles = np.random.randint(1, n_haplotypes // 2)
        alt_haplotype_indices = np.random.choice(range(n_haplotypes), num_alt_alleles, replace=False)
        true_haplotypes[alt_haplotype_indices, site] = 1

    # --- 3. Calculate True Genetic Diversity Metrics ---
    # a_n is the correction factor for Watterson's estimator
    a_n = sum(1.0 / i for i in range(1, n_haplotypes))

    # S is the number of segregating sites
    S_true = np.sum(np.sum(true_haplotypes, axis=0) > 0)

    # Watterson's estimator (theta)
    theta_w_true = S_true / a_n

    # Nucleotide diversity (pi)
    p_true = np.mean(true_haplotypes, axis=0) # frequency of alt allele at each site
    q_true = 1 - p_true
    # pi is the sum of heterozygosity (2*p*q) over all sites
    pi_true = np.sum(2 * p_true * q_true)

    print("--- 1. True Population Values (Before Filtering) ---")
    print(f"True number of segregating sites (S): {S_true}")
    print(f"True Watterson's estimator (theta): {theta_w_true:.4f}")
    print(f"True Nucleotide Diversity (pi): {pi_true:.4f}\n")


    # --- 4. Simulate Filtering and Imputation ---
    observed_haplotypes = true_haplotypes.copy()
    # Find all alternative alleles (value of 1)
    alt_allele_locations = np.where(observed_haplotypes == 1)
    
    # Iterate through each alt allele and decide if it gets filtered
    num_filtered = 0
    for hap_idx, site_idx in zip(*alt_allele_locations):
        if np.random.random() < p_filter:
            # Filter and impute with reference allele (0)
            observed_haplotypes[hap_idx, site_idx] = 0
            num_filtered += 1
            
    print(f"--- 2. Simulating the Lab Procedure ---")
    print(f"A total of {num_filtered} alternative alleles were filtered and replaced with the reference allele.\n")

    # --- 5. Calculate Observed Genetic Diversity Metrics ---
    S_observed = np.sum(np.sum(observed_haplotypes, axis=0) > 0)
    theta_w_observed = S_observed / a_n
    
    p_observed = np.mean(observed_haplotypes, axis=0)
    q_observed = 1 - p_observed
    pi_observed = np.sum(2 * p_observed * q_observed)

    print("--- 3. Observed Values (After Filtering and Imputation) ---")
    print(f"Observed number of segregating sites (S): {S_observed}")
    print(f"Observed Watterson's estimator (theta): {theta_w_observed:.4f}")
    print(f"Observed Nucleotide Diversity (pi): {pi_observed:.4f}\n")
    
    # --- 6. Conclusion ---
    print("--- Conclusion ---")
    print("The filtering and imputation process causes some true polymorphic sites to appear monomorphic,")
    print(f"reducing S from {S_true} to {S_observed}. This systematically lowers theta.")
    print("\nThis process also reduces the frequency of the alternative allele at the remaining polymorphic sites.")
    print(f"This reduces overall heterozygosity, causing pi to be underestimated (from {pi_true:.4f} to {pi_observed:.4f}).")
    print("\nTherefore, both Watterson's estimator (theta) and nucleotide diversity (pi) are biased.")

if __name__ == '__main__':
    run_simulation()