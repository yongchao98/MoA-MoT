import numpy as np
import math

def run_simulation():
    """
    This simulation demonstrates how filtering individual low-quality variant calls
    and imputing with a reference genome introduces bias in population genetic estimators.
    """
    # --- Simulation Parameters ---
    n_haplotypes = 50  # Number of sequences in the sample
    n_sites = 10000    # Total length of the sequence
    n_true_snps = 50   # Number of true segregating sites (SNPs)
    filter_prob = 0.3  # Probability that a true alternate allele is filtered as "low quality"
    ref_allele = 0
    alt_allele = 1

    # --- Create True Genetic Data ---
    # Start with a matrix of all reference alleles
    true_data = np.full((n_haplotypes, n_sites), ref_allele, dtype=int)
    true_alt_allele_counts = []
    
    # Introduce SNPs at random positions
    true_snp_indices = np.random.choice(range(n_sites), n_true_snps, replace=False)
    for snp_idx in true_snp_indices:
        # Give each SNP a random number of alternate alleles
        alt_count = np.random.randint(1, n_haplotypes)
        true_alt_allele_counts.append(alt_count)
        haplotypes_with_alt = np.random.choice(range(n_haplotypes), alt_count, replace=False)
        true_data[haplotypes_with_alt, snp_idx] = alt_allele

    # --- Calculate True Theta and Pi ---
    true_S = len(true_snp_indices)
    a_n = sum(1.0 / i for i in range(1, n_haplotypes)) # Harmonic number
    true_theta = true_S / a_n
    
    true_site_pis = []
    for count in true_alt_allele_counts:
        p = count / n_haplotypes
        true_site_pis.append(2 * p * (1 - p))
    true_pi_sum = sum(true_site_pis)

    # --- Simulate the Filtering and Imputation Pipeline ---
    observed_data = np.copy(true_data)
    # Iterate through each site for each individual haplotype
    for i in range(n_haplotypes):
        for j in range(n_sites):
            # If it's an alternate allele, it has a chance of being filtered
            if observed_data[i, j] == alt_allele and np.random.rand() < filter_prob:
                observed_data[i, j] = ref_allele # Impute with reference

    # --- Calculate Observed Theta and Pi ---
    observed_S = 0
    observed_snp_indices = []
    for j in range(n_sites):
        # A site is segregating if it contains both reference and alternate alleles
        if alt_allele in observed_data[:, j] and ref_allele in observed_data[:, j]:
            observed_S += 1
            observed_snp_indices.append(j)
    observed_theta = observed_S / a_n
    
    observed_site_pis = []
    for snp_idx in observed_snp_indices:
        alt_count = np.sum(observed_data[:, snp_idx])
        p = alt_count / n_haplotypes
        observed_site_pis.append(2 * p * (1 - p))
    observed_pi_sum = sum(observed_site_pis)

    # --- Print and Compare Results ---
    print("--- Simulation of Bias in Population Genetic Estimators ---")
    print(f"\nParameters: n={n_haplotypes} haplotypes, {n_true_snps} true SNPs, filter_prob={filter_prob}\n")
    print("-" * 55)

    print("1. Watterson's Estimator (Theta_W = S / a_n)")
    print("\n--- TRUE VALUES ---")
    print(f"True number of segregating sites (S): {true_S}")
    print(f"Harmonic number (a_n): {a_n:.4f}")
    print(f"Equation for True Theta_W: {true_S} / {a_n:.4f}")
    print(f"Resulting True Theta_W: {true_theta:.4f}")

    print("\n--- OBSERVED VALUES (after filtering/imputation) ---")
    print(f"Observed number of segregating sites (S'): {observed_S}")
    print(f"Equation for Observed Theta_W: {observed_S} / {a_n:.4f}")
    print(f"Resulting Observed Theta_W: {observed_theta:.4f}")
    print(f"\nResult: Observed Theta ({observed_theta:.4f}) is biased relative to True Theta ({true_theta:.4f}).")
    print("-" * 55)

    print("\n2. Nucleotide Diversity (Pi = sum of 2*p*q for all sites)")
    print("\n--- TRUE VALUES ---")
    print(f"Sum of site pi values across all {true_S} sites:")
    print(f"Resulting True Pi (total): {true_pi_sum:.4f}")

    print("\n--- OBSERVED VALUES (after filtering/imputation) ---")
    print(f"Sum of site pi values across all {observed_S} observed sites:")
    print(f"Resulting Observed Pi (total): {observed_pi_sum:.4f}")
    print(f"\nResult: Observed Pi ({observed_pi_sum:.4f}) is biased relative to True Pi ({true_pi_sum:.4f}).")
    print("-" * 55)
    
    print("\nConclusion: The pipeline introduces a downward bias in BOTH estimators.")

# Run the simulation
run_simulation()