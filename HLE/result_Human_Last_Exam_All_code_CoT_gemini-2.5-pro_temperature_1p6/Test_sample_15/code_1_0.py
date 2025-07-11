import msprime
import numpy as np

def calculate_pi_from_freq(genotype_matrix):
    """
    Calculates nucleotide diversity (pi) from a genotype matrix using allele frequencies.
    This is an unbiased estimator for the population diversity.
    """
    n_haplotypes, n_sites = genotype_matrix.shape
    if n_haplotypes < 2:
        return 0.0
        
    # Calculate pi for each site: 2 * p * q where p and q are allele frequencies
    # Sum across all sites to get total pi for the sequence.
    alt_allele_counts = np.sum(genotype_matrix, axis=0)
    alt_allele_freqs = alt_allele_counts / n_haplotypes
    
    # pi is the sum of heterozygosity over all sites
    pi_sum = np.sum(2 * alt_allele_freqs * (1 - alt_allele_freqs))
    
    # Correction for sample size to make it an unbiased estimator
    pi = pi_sum * n_haplotypes / (n_haplotypes - 1)
    
    return pi

def run_simulation():
    """
    Simulates the bioinformatics pipeline and demonstrates the bias in
    Watterson's theta and nucleotide diversity pi.
    """
    # --- 1. Simulation Parameters ---
    num_samples = 50      # Number of diploid individuals
    num_haplotypes = 2 * num_samples
    seq_length = 500_000  # Length of the sequence
    random_seed = 123     # For reproducibility
    # Fraction of a sample's variants to filter out due to "low quality"
    filter_fraction = 0.25 # 25%

    print("Step 1: Simulating 'true' genetic data for a population...")
    # Generate a true genetic history and mutations using msprime
    ts_true = msprime.sim_ancestry(
        samples=num_samples,
        population_size=10_000,
        sequence_length=seq_length,
        random_seed=random_seed
    )
    ts_true = msprime.sim_mutations(ts_true, rate=1e-8, random_seed=random_seed)

    # --- 2. Calculate True Statistics ---
    S_true = ts_true.num_sites
    pi_true = ts_true.diversity()
    a_n = sum(1.0 / i for i in range(1, num_haplotypes))
    theta_w_true = S_true / a_n

    print("\n--- Ground Truth Statistics ---")
    print(f"Number of haplotypes (2n): {num_haplotypes}")
    print(f"True number of segregating sites (S_true): {S_true}")
    print(f"True Watterson's estimator (theta_w_true) based on S_true: ")
    print(f"  theta_w = {S_true} / {a_n:.2f} = {theta_w_true:.4f}")
    print(f"True nucleotide diversity (pi_true): {pi_true:.4f}\n")

    # --- 3. Simulate the Flawed Filtering and Imputation Process ---
    print("Step 2: Simulating the filtering and imputation process...")
    print(f"For each of the {num_samples} samples, {filter_fraction*100}% of its variant sites will be removed and imputed with the reference allele.")
    
    # Get the genotype matrix (haplotypes x sites)
    genotype_matrix_true = ts_true.genotype_matrix()
    genotype_matrix_observed = np.copy(genotype_matrix_true)

    for i in range(num_samples):
        # Define the two haplotypes for this diploid sample
        h1_idx, h2_idx = 2 * i, 2 * i + 1
        
        # Find sites where this specific sample has at least one variant allele
        sample_snv_indices = np.where(
            (genotype_matrix_true[h1_idx, :] == 1) | (genotype_matrix_true[h2_idx, :] == 1)
        )[0]
        
        num_to_filter = int(len(sample_snv_indices) * filter_fraction)

        if num_to_filter > 0:
            # Randomly choose sites to filter for this sample
            indices_to_filter = np.random.choice(
                sample_snv_indices, size=num_to_filter, replace=False
            )
            # Impute with reference (0) for both haplotypes at these sites
            genotype_matrix_observed[h1_idx, indices_to_filter] = 0
            genotype_matrix_observed[h2_idx, indices_to_filter] = 0

    # --- 4. Calculate Observed Statistics ---
    print("\nStep 3: Calculating statistics from the processed data...")

    # Observed S is the number of sites that are still polymorphic
    site_sums = np.sum(genotype_matrix_observed, axis=0)
    S_observed = np.count_nonzero(site_sums)

    # Observed pi is calculated from the modified matrix
    pi_observed = calculate_pi_from_freq(genotype_matrix_observed)

    # Observed Watterson's theta
    theta_w_observed = S_observed / a_n

    print("\n--- Observed Statistics (After Flawed Processing) ---")
    print(f"Observed number of segregating sites (S_observed): {S_observed}")
    print(f"Observed Watterson's estimator (theta_w_observed):")
    print(f"  theta_w = {S_observed} / {a_n:.2f} = {theta_w_observed:.4f}")
    print(f"Observed nucleotide diversity (pi_observed): {pi_observed:.4f}\n")
    
    # --- 5. Final Analysis and Conclusion ---
    print("--- Analysis of Bias ---")
    print("Bias in Watterson's Estimator (theta):")
    print("The imputation process converts true variant alleles to reference alleles. When all variant alleles at a site are removed across the entire sample (e.g., if it was a rare variant), that site no longer appears to be segregating.")
    print(f"This leads to an underestimation of the number of segregating sites (S_observed = {S_observed} vs. S_true = {S_true}).")
    print("Since S is underestimated, Watterson's estimator is biased downwards.\n")

    print("Bias in Nucleotide Diversity (pi):")
    print("Pi measures the average genetic difference between sequences. By replacing true variants with the reference allele, sequences are made artificially more similar to each other.")
    print("This happens because the frequency of variant alleles is lowered, which reduces the average number of pairwise differences.")
    print(f"This leads to an underestimation of nucleotide diversity (pi_observed = {pi_observed:.4f} vs. pi_true = {pi_true:.4f}).")
    print("Therefore, pi is also biased downwards.\n")
    
    print("Conclusion: The described workflow biases BOTH Watterson's estimator and nucleotide diversity.")

if __name__ == '__main__':
    run_simulation()