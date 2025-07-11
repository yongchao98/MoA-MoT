import numpy as np
import random

def calculate_s(data):
    """Calculates the number of segregating sites (S)."""
    n_haplotypes, n_sites = data.shape
    s = 0
    for j in range(n_sites):
        site_column = data[:, j]
        # A site is segregating if it has both reference (0) and alternative (1) alleles
        if 0 < np.sum(site_column) < n_haplotypes:
            s += 1
    return s

def calculate_pi_components(data):
    """Calculates components for pi: total pairwise differences and number of pairs."""
    n_haplotypes, n_sites = data.shape
    if n_haplotypes < 2:
        return 0, 1
    
    total_pairwise_diffs = 0
    # The number of unique pairs of sequences is n*(n-1)/2
    num_pairs = n_haplotypes * (n_haplotypes - 1) / 2
    
    for j in range(n_sites):
        # Count the alternative allele (1) at this site
        k = np.sum(data[:, j])
        # The number of pairs that differ at this site is k * (n_haplotypes - k)
        diffs_at_site = k * (n_haplotypes - k)
        total_pairwise_diffs += diffs_at_site
        
    return total_pairwise_diffs, num_pairs

# --- Simulation Parameters ---
n_haplotypes = 50
n_sites = 10000
n_polymorphic_sites = 200
# Probability that a true alternative allele call is filtered and imputed with reference
filter_prob = 0.30

# --- 1. Generate True Data ---
print("--- Step 1: Simulating True Genetic Data ---")
# Start with all reference alleles (0)
true_data = np.zeros((n_haplotypes, n_sites), dtype=int)

# Select sites to be polymorphic
polymorphic_indices = random.sample(range(n_sites), n_polymorphic_sites)

for site_idx in polymorphic_indices:
    # For each polymorphic site, decide how many haplotypes carry the alternative allele
    alt_allele_count = random.randint(1, n_haplotypes - 1)
    # Randomly assign the alternative alleles (1) to haplotypes
    haplotype_indices_with_alt = random.sample(range(n_haplotypes), alt_allele_count)
    true_data[haplotype_indices_with_alt, site_idx] = 1
print(f"Generated a dataset with {n_haplotypes} samples and {n_polymorphic_sites} true polymorphic sites.\n")

# --- 2. Calculate True Statistics ---
print("--- Step 2: Calculating Statistics from True Data ---")
# Calculate harmonic number a_n = sum_{i=1}^{n-1} 1/i
a_n = np.sum(1.0 / np.arange(1, n_haplotypes))

# Watterson's Theta
true_S = calculate_s(true_data)
true_theta = true_S / a_n
print("Watterson's Theta (based on true data):")
print(f"Formula: S / a_n")
print(f"Result: {true_S} / {a_n:.4f} = {true_theta:.4f}\n")

# Nucleotide Diversity (pi)
total_diffs_true, num_pairs = calculate_pi_components(true_data)
true_pi = total_diffs_true / num_pairs
print("Nucleotide Diversity (pi) (based on true data):")
print(f"Formula: Total_Pairwise_Differences / Number_of_Pairs")
print(f"Result: {total_diffs_true} / {int(num_pairs)} = {true_pi:.4f}\n")


# --- 3. Simulate Filtering and Imputation ---
print(f"--- Step 3: Simulating Filtering ({filter_prob*100}% probability) and Imputation ---")
observed_data = np.copy(true_data)
# Iterate through the data and apply filtering/imputation
for i in range(n_haplotypes):
    for j in range(n_sites):
        # If the true allele is the alternative allele (1)...
        if true_data[i, j] == 1:
            # ...there's a chance it gets filtered and imputed with the reference (0)
            if random.random() < filter_prob:
                observed_data[i, j] = 0
print("Pipeline simulation complete.\n")

# --- 4. Calculate Observed (Biased) Statistics ---
print("--- Step 4: Calculating Statistics from Observed (Post-Pipeline) Data ---")
# Watterson's Theta
observed_S = calculate_s(observed_data)
observed_theta = observed_S / a_n
print("Watterson's Theta (based on observed data):")
print(f"Formula: S / a_n")
print(f"Result: {observed_S} / {a_n:.4f} = {observed_theta:.4f}\n")

# Nucleotide Diversity (pi)
total_diffs_obs, num_pairs = calculate_pi_components(observed_data)
observed_pi = total_diffs_obs / num_pairs
print("Nucleotide Diversity (pi) (based on observed data):")
print(f"Formula: Total_Pairwise_Differences / Number_of_Pairs")
print(f"Result: {total_diffs_obs} / {int(num_pairs)} = {observed_pi:.4f}\n")

# --- 5. Conclusion ---
print("--- Conclusion ---")
print(f"The true number of polymorphic sites was {true_S}, but after the pipeline, we only observe {observed_S}.")
print(f"The true theta was {true_theta:.4f}, but the observed theta is {observed_theta:.4f}.")
print(f"The true pi was {true_pi:.4f}, but the observed pi is {observed_pi:.4f}.")
print("Both estimators are biased downwards due to the systematic removal of alternative alleles.")
