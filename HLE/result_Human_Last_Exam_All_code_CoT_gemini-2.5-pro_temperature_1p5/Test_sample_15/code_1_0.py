import numpy as np
import random

def calculate_pi(haplotypes):
    """Calculates nucleotide diversity (pi) as the mean of 2*p*q across sites."""
    n_haplotypes, n_sites = haplotypes.shape
    if n_haplotypes < 2:
        return 0.0
    
    # Calculate allele frequencies at each site
    alt_counts = np.sum(haplotypes, axis=0) # Sum of alternate alleles (1s) per site
    p_alt = alt_counts / n_haplotypes
    p_ref = 1.0 - p_alt
    
    # Pi at each site is 2 * p_ref * p_alt
    pi_per_site = 2 * p_ref * p_alt
    
    # Overall pi is the average across all sites
    return np.mean(pi_per_site)

def calculate_S(haplotypes):
    """Calculates the number of segregating sites (S)."""
    n_haplotypes = haplotypes.shape[0]
    alt_counts = np.sum(haplotypes, axis=0)
    # A site is segregating if the alternate allele count is between 1 and n-1
    segregating_sites = np.where((alt_counts > 0) & (alt_counts < n_haplotypes))[0]
    return len(segregating_sites)
    
def calculate_watterson_theta(S, n_haplotypes):
    """Calculates Watterson's theta estimator (per site)."""
    if n_haplotypes < 2:
        return 0.0
    # a_n is the (n-1)th harmonic number
    a_n = np.sum(1.0 / np.arange(1, n_haplotypes))
    return S / a_n

# --- Simulation Parameters ---
n_samples = 50       # Number of diploid individuals
n_haplotypes = n_samples * 2
n_sites = 20000      # Length of the sequence
true_n_snps = 300    # Number of true SNVs
p_filter = 0.25      # 25% chance a true alternate allele call is filtered

# --- 1. Simulate True Genetic Data ---
# Start with all reference alleles (0)
true_haplotypes = np.zeros((n_haplotypes, n_sites), dtype=int)

# Randomly select sites to be SNVs
snp_positions = sorted(random.sample(range(n_sites), true_n_snps))

# For each SNV, distribute alternate alleles (1) to create a realistic
# site frequency spectrum (SFS) with more rare variants.
print("Step 1: Simulating true genetic data...")
for pos in snp_positions:
    # Use a beta distribution to favor low-frequency variants
    allele_freq = np.random.beta(a=0.5, b=10)
    n_alt_alleles = int(round(allele_freq * n_haplotypes))
    if n_alt_alleles == 0: 
        n_alt_alleles = 1 # Ensure site is truly segregating
    
    alt_haplotype_indices = random.sample(range(n_haplotypes), n_alt_alleles)
    true_haplotypes[alt_haplotype_indices, pos] = 1

# --- 2. Calculate Statistics from True Data ---
print("\nStep 2: Calculating statistics from the original, correct data...")
S_true = calculate_S(true_haplotypes)
pi_true = calculate_pi(true_haplotypes)
# We calculate theta per site to compare it with pi
theta_true_raw = calculate_watterson_theta(S_true, n_haplotypes)
theta_true_per_site = theta_true_raw / n_sites

print(f"--- True Data ---")
print(f"Number of Segregating Sites (S_true): {S_true}")
print(f"Nucleotide Diversity (pi_true): {pi_true:.6f}")
print(f"Watterson's Estimator per site (theta_true): {theta_true_per_site:.6f}")


# --- 3. Simulate Flawed Processing Pipeline ---
print(f"\nStep 3: Simulating the flawed pipeline (filtering + imputation)...")
observed_haplotypes = np.copy(true_haplotypes)
# Find every instance of an alternate allele
alt_allele_indices = np.where(observed_haplotypes == 1)

# Iterate through each true alternate allele and decide if it gets filtered
filtered_count = 0
for i in range(len(alt_allele_indices[0])):
    if random.random() < p_filter:
        haplotype_idx = alt_allele_indices[0][i]
        site_idx = alt_allele_indices[1][i]
        # Impute with reference allele
        observed_haplotypes[haplotype_idx, site_idx] = 0
        filtered_count += 1
print(f"Filtered and replaced {filtered_count} alternate alleles with reference alleles.")

# --- 4. Calculate Statistics from Observed Data ---
print("\nStep 4: Calculating statistics from the processed (flawed) data...")
S_obs = calculate_S(observed_haplotypes)
pi_obs = calculate_pi(observed_haplotypes)
theta_obs_raw = calculate_watterson_theta(S_obs, n_haplotypes)
theta_obs_per_site = theta_obs_raw / n_sites

print(f"--- Observed Data ---")
print(f"Number of Segregating Sites (S_obs): {S_obs}")
print(f"Nucleotide Diversity (pi_obs): {pi_obs:.6f}")
print(f"Watterson's Estimator per site (theta_obs): {theta_obs_per_site:.6f}")

# --- 5. Final Conclusion ---
print("\n--- Final Conclusion ---")
s_bias = (S_obs - S_true) / S_true * 100 if S_true > 0 else 0
pi_bias = (pi_obs - pi_true) / pi_true * 100 if pi_true > 0 else 0

print(f"Bias in S: The number of segregating sites is underestimated. (Observed: {S_obs}, True: {S_true}, Bias: {s_bias:.2f}%)")
print("--> This causes Watterson's estimator (theta) to be BIASED.\n")

print(f"Bias in pi: Nucleotide diversity is underestimated. (Observed: {pi_obs:.6f}, True: {pi_true:.6f}, Bias: {pi_bias:.2f}%)")
print("--> This causes nucleotide diversity (pi) to be BIASED.\n")

print("Both Watterson's estimator (theta) and pi (nucleotide diversity) are biased.")
