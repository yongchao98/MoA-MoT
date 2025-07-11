import numpy as np
import math

def calculate_s(haplotypes):
    """Calculates the number of segregating (polymorphic) sites."""
    n_sites = haplotypes.shape[1]
    s = 0
    for i in range(n_sites):
        # A site is segregating if more than one allele type exists
        if len(np.unique(haplotypes[:, i])) > 1:
            s += 1
    return s

def calculate_pi(haplotypes):
    """Calculates nucleotide diversity (pi) as the average pairwise difference."""
    n_haplotypes, n_sites = haplotypes.shape
    if n_haplotypes < 2:
        return 0.0
    
    total_pairwise_diffs = 0
    num_pairs = 0
    
    # Iterate through all unique pairs of haplotypes
    for i in range(n_haplotypes):
        for j in range(i + 1, n_haplotypes):
            # Sum of differences between haplotype i and haplotype j
            diffs = np.sum(haplotypes[i, :] != haplotypes[j, :])
            total_pairwise_diffs += diffs
            num_pairs += 1
            
    # Average number of differences per site
    pi = total_pairwise_diffs / num_pairs
    return pi

def run_simulation():
    """
    Simulates the bioinformatics pipeline and calculates the bias in theta and pi.
    """
    # --- 1. Simulation Parameters ---
    n_haplotypes = 100  # e.g., 50 diploid samples
    n_sites = 20000
    n_true_polymorphic_sites = 100
    # Probability that a true alternate allele is filtered and imputed as reference
    filter_and_impute_prob = 0.25 

    # Allele representation: 0 = reference, 1 = alternate
    ref_allele = 0
    alt_allele = 1

    print("--- Simulation Setup ---")
    print(f"Number of haplotypes: {n_haplotypes}")
    print(f"Number of sites: {n_sites}")
    print(f"True number of polymorphic sites: {n_true_polymorphic_sites}")
    print(f"Probability of filtering an alternate allele: {filter_and_impute_prob}\n")

    # --- 2. Generate 'True' Haplotype Data ---
    # Start with a matrix of all reference alleles
    true_haplotypes = np.full((n_haplotypes, n_sites), ref_allele, dtype=int)

    # Introduce polymorphic sites
    polymorphic_indices = np.random.choice(n_sites, n_true_polymorphic_sites, replace=False)
    for site_idx in polymorphic_indices:
        # For each polymorphic site, give it a random minor allele frequency (MAF)
        maf = np.random.uniform(0.01, 0.5)
        n_alt_alleles = int(round(maf * n_haplotypes))
        alt_haplotype_indices = np.random.choice(n_haplotypes, n_alt_alleles, replace=False)
        true_haplotypes[alt_haplotype_indices, site_idx] = alt_allele

    # --- 3. Calculate True Genetic Diversity Metrics ---
    s_true = calculate_s(true_haplotypes)
    pi_true = calculate_pi(true_haplotypes)
    # Watterson's estimator constant
    a_n = sum(1.0 / i for i in range(1, n_haplotypes))
    theta_w_true = s_true / a_n

    print("--- 'True' Genetic Diversity (Before Filtering) ---")
    print(f"True Number of Segregating Sites (S): {s_true}")
    print(f"True Nucleotide Diversity (pi): {pi_true:.4f}")
    print(f"True Watterson's Estimator (theta): {theta_w_true:.4f}\n")

    # --- 4. Simulate Filtering and Imputation ---
    observed_haplotypes = np.copy(true_haplotypes)
    # Find all occurrences of the alternate allele
    alt_allele_locations = np.where(observed_haplotypes == alt_allele)
    
    for i in range(len(alt_allele_locations[0])):
        # For each alternate allele, there's a chance it gets filtered
        if np.random.rand() < filter_and_impute_prob:
            row = alt_allele_locations[0][i]
            col = alt_allele_locations[1][i]
            # Impute with the reference allele
            observed_haplotypes[row, col] = ref_allele

    # --- 5. Calculate Observed Genetic Diversity Metrics ---
    s_observed = calculate_s(observed_haplotypes)
    pi_observed = calculate_pi(observed_haplotypes)
    theta_w_observed = s_observed / a_n

    print("--- 'Observed' Genetic Diversity (After Filtering & Imputation) ---")
    print(f"Observed Number of Segregating Sites (S): {s_observed}")
    print(f"Observed Nucleotide Diversity (pi): {pi_observed:.4f}")
    print(f"Observed Watterson's Estimator (theta): {theta_w_observed:.4f}\n")

    # --- 6. Conclusion ---
    print("--- Conclusion ---")
    s_bias = "biased downwards" if s_observed < s_true else "not biased"
    pi_bias = "biased downwards" if pi_observed < pi_true else "not biased"
    theta_bias = "biased downwards" if theta_w_observed < theta_w_true else "not biased"

    print(f"Comparison shows S is {s_bias} ({s_observed} < {s_true}).")
    print(f"This causes Watterson's theta to be {theta_bias} ({theta_w_observed:.4f} < {theta_w_true:.4f}).")
    print(f"Comparison shows pi is {pi_bias} ({pi_observed:.4f} < {pi_true:.4f}).")
    print("\nTherefore, both Watterson's estimator (theta) and pi (nucleotide diversity) are biased.")

if __name__ == '__main__':
    run_simulation()