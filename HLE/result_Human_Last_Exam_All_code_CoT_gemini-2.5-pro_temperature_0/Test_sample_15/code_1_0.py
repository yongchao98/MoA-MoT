import numpy as np
import random

def run_simulation():
    """
    Simulates the bioinformatics pipeline to demonstrate its effect on
    Watterson's theta and nucleotide diversity (pi).
    """
    # --- 1. Parameters ---
    num_individuals = 50
    num_haplotypes = num_individuals * 2
    seq_len = 20000
    true_num_snps = 100  # The true number of segregating sites (S)
    filter_prob = 0.3  # 30% chance a variant call is low quality and gets filtered

    # --- Helper Functions ---
    def calculate_pi(haplotypes):
        n = haplotypes.shape[0]
        if n < 2: return 0.0
        num_pairwise_comparisons = n * (n - 1) / 2
        num_alt = np.sum(haplotypes, axis=0)
        num_ref = n - num_alt
        total_pairwise_diffs = np.sum(num_ref * num_alt)
        return total_pairwise_diffs, num_pairwise_comparisons

    def calculate_s(haplotypes):
        n = haplotypes.shape[0]
        column_sums = np.sum(haplotypes, axis=0)
        return np.sum((column_sums > 0) & (column_sums < n))

    # --- 2. Generate True Data ---
    true_haplotypes = np.zeros((num_haplotypes, seq_len), dtype=int)
    snp_positions = random.sample(range(seq_len), true_num_snps)
    for pos in snp_positions:
        maf = random.uniform(0.01, 0.5) # Use a range of allele frequencies
        num_variants = max(1, int(round(maf * num_haplotypes)))
        variant_indices = random.sample(range(num_haplotypes), num_variants)
        true_haplotypes[variant_indices, pos] = 1

    # --- 3. Calculate True Statistics ---
    true_S = calculate_s(true_haplotypes)
    a_n = sum(1.0 / i for i in range(1, num_haplotypes))
    true_theta = true_S / a_n
    true_total_diffs, num_comparisons = calculate_pi(true_haplotypes)
    true_pi = true_total_diffs / num_comparisons

    # --- 4. Simulate Lab Pipeline ---
    observed_haplotypes = np.copy(true_haplotypes)
    variant_locations = np.where(observed_haplotypes == 1)
    for hap_idx, site_idx in zip(*variant_locations):
        if random.random() < filter_prob:
            observed_haplotypes[hap_idx, site_idx] = 0 # Impute with reference

    # --- 5. Calculate Observed Statistics ---
    observed_S = calculate_s(observed_haplotypes)
    observed_theta = observed_S / a_n
    observed_total_diffs, _ = calculate_pi(observed_haplotypes)
    observed_pi = observed_total_diffs / num_comparisons

    # --- 6. Print Results and Explanation ---
    print("This simulation demonstrates the bias introduced by the described pipeline.")
    
    # --- Watterson's Theta Analysis ---
    print("\n--- Analysis of Watterson's Estimator (theta) ---")
    print("Equation: theta = S / a_n")
    print(f"Based on n = {num_haplotypes} haplotypes, the constant a_n is {a_n:.4f}")
    print("\n[TRUE CALCULATION]")
    print(f"True number of segregating sites (S_true): {true_S}")
    print(f"theta_true = {true_S} / {a_n:.4f} = {true_theta:.4f}")
    print("\n[OBSERVED CALCULATION]")
    print(f"Observed number of segregating sites (S_obs): {observed_S}")
    print(f"theta_obs = {observed_S} / {a_n:.4f} = {observed_theta:.4f}")
    print("Conclusion: Theta is biased because the filtering of rare variants can cause S to be underestimated.")

    # --- Pi Analysis ---
    print("\n--- Analysis of Nucleotide Diversity (pi) ---")
    print("Equation: pi = (Total Pairwise Differences) / (Number of Pairs)")
    print(f"Number of haplotype pairs is {int(num_comparisons)}")
    print("\n[TRUE CALCULATION]")
    print(f"True total pairwise differences: {int(true_total_diffs)}")
    print(f"pi_true = {int(true_total_diffs)} / {int(num_comparisons)} = {true_pi:.4f}")
    print("\n[OBSERVED CALCULATION]")
    print(f"Observed total pairwise differences: {int(observed_total_diffs)}")
    print(f"pi_obs = {int(observed_total_diffs)} / {int(num_comparisons)} = {observed_pi:.4f}")
    print("Conclusion: Pi is biased because imputing filtered variants with the reference allele systematically lowers variant frequencies.")

    print("\n--- OVERALL CONCLUSION ---")
    print("Both Watterson's estimator (theta) and pi (nucleotide diversity) are biased by this pipeline.")

if __name__ == '__main__':
    run_simulation()