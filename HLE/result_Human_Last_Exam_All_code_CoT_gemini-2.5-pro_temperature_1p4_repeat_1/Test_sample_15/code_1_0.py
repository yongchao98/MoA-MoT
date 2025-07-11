import numpy as np
import random

def run_simulation_and_calculate_bias():
    """
    Simulates a bioinformatics pipeline to demonstrate bias in population genetics estimators.

    This function performs the following steps:
    1. Defines parameters for the simulation (e.g., number of haplotypes, sequence length).
    2. Generates a 'true' dataset of haplotypes with a known number of segregating sites.
    3. Simulates a filtering/imputation process where true variants are randomly replaced
       with the reference allele, creating an 'observed' dataset.
    4. Calculates Watterson's theta and nucleotide diversity (pi) for both datasets.
    5. Prints the results to show the systematic underestimation (bias) in the observed values.
    """
    # 1. Define simulation parameters
    num_haplotypes = 100  # n, number of sequences in the sample
    seq_length = 50000    # L, total length of the sequence region
    true_s = 200          # S, the true number of segregating (polymorphic) sites
    filter_prob = 0.30    # Probability of a variant call being filtered for a given sample

    print(f"--- Simulation Parameters ---")
    print(f"Number of Haplotypes (n): {num_haplotypes}")
    print(f"True Number of Segregating Sites (S): {true_s}")
    print(f"Filtering Probability per Variant: {filter_prob}\n")

    # 2. Generate 'True' Haplotype Data
    # Start with a matrix of all reference alleles (0)
    true_haplotypes = np.zeros((num_haplotypes, seq_length), dtype=int)
    
    # Select sites to be polymorphic
    polymorphic_sites = random.sample(range(seq_length), true_s)
    
    # For each polymorphic site, introduce variant alleles (1)
    for site in polymorphic_sites:
        # For simplicity, assign a random number of variant alleles (from 1 to n-1)
        num_variants = random.randint(1, num_haplotypes - 1)
        variant_haplotypes = random.sample(range(num_haplotypes), num_variants)
        true_haplotypes[variant_haplotypes, site] = 1

    # 3. Simulate Filtering and Imputation to create 'Observed' data
    observed_haplotypes = np.copy(true_haplotypes)
    variants_filtered = 0
    
    # Iterate through each true variant and apply the filter
    for i in range(num_haplotypes):
        for j in range(seq_length):
            if observed_haplotypes[i, j] == 1:
                if random.random() < filter_prob:
                    # Filtered variant is imputed with the reference allele (0)
                    observed_haplotypes[i, j] = 0
                    variants_filtered += 1
    
    print(f"--- Pipeline Simulation ---")
    print(f"A total of {variants_filtered} individual variant alleles were filtered out and replaced with the reference allele.\n")

    # 4. Calculate Estimators for both 'True' and 'Observed' data
    
    # --- Helper function to calculate Pi (average pairwise differences) ---
    def calculate_pi(haplotypes):
        n = haplotypes.shape[0]
        num_pairs = n * (n - 1) / 2
        total_diffs = 0
        # Iterate over all sites
        for j in range(haplotypes.shape[1]):
            # Count number of variant alleles (k) at the site
            k = np.sum(haplotypes[:, j])
            # Number of differences at this site = k * (n - k)
            total_diffs += k * (n - k)
        return total_diffs / num_pairs

    # --- Helper function to calculate S (number of segregating sites) ---
    def calculate_s(haplotypes):
        # A site is segregating if the sum of alleles in its column is > 0 and < n
        # Since we only have 0s and 1s, we just need to check if the sum > 0
        site_sums = np.sum(haplotypes, axis=0)
        return np.sum(site_sums > 0)

    # --- Calculations ---
    # S
    s_true = calculate_s(true_haplotypes)
    s_observed = calculate_s(observed_haplotypes)

    # Pi
    pi_true = calculate_pi(true_haplotypes)
    pi_observed = calculate_pi(observed_haplotypes)

    # Watterson's Theta
    # First, calculate the scaling factor a_n
    a_n = sum(1.0 / i for i in range(1, num_haplotypes))
    theta_w_true = s_true / a_n
    theta_w_observed = s_observed / a_n

    # 5. Print Results
    print("--- Estimator Calculations ---")
    print("\n--- Nucleotide Diversity (Pi) ---")
    print(f"Pi (True):     {pi_true:.4f} (Calculated from original data)")
    print(f"Pi (Observed): {pi_observed:.4f} (Calculated after filtering/imputation)")
    if pi_observed < pi_true:
        print("Result: Observed Pi is lower than True Pi. The estimator is BIASED.")
    else:
        print("Result: No bias detected in this run (statistically unlikely).")

    print("\n--- Watterson's Estimator (Theta) ---")
    print(f"Number of Segregating Sites (S_true):     {s_true}")
    print(f"Number of Segregating Sites (S_observed): {s_observed}")
    print(f"Harmonic number (a_n for n={num_haplotypes}): {a_n:.4f}")
    
    print("\nFinal Equation for Theta (True):")
    print(f"Theta_W = S / a_n = {s_true} / {a_n:.4f} = {theta_w_true:.4f}")

    print("\nFinal Equation for Theta (Observed):")
    print(f"Theta_W = S / a_n = {s_observed} / {a_n:.4f} = {theta_w_observed:.4f}")

    if theta_w_observed < theta_w_true:
        print("\nResult: Observed Theta is lower than True Theta. The estimator is BIASED.")
    else:
        print("\nResult: No bias detected in this run (statistically unlikely).")

    print("\n--- Overall Conclusion ---")
    print("The simulation shows that the data processing pipeline systematically reduces both the number of segregating sites (S) and the average pairwise differences (Pi).")
    print("Therefore, both Watterson's estimator (theta) and nucleotide diversity (pi) are biased.")

# Run the simulation
run_simulation_and_calculate_bias()
<<<C>>>