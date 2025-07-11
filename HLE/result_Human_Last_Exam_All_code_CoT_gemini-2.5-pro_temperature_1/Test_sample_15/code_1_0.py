import random
import numpy as np

def calculate_pi(haplotypes):
    """Calculates nucleotide diversity (pi) for a set of haplotypes."""
    n_haplotypes = len(haplotypes)
    if n_haplotypes < 2:
        return 0.0
    
    seq_len = len(haplotypes[0])
    pairwise_diffs_sum = 0
    num_pairs = n_haplotypes * (n_haplotypes - 1) / 2

    # This is more efficient than direct pairwise comparison for long sequences
    for i in range(seq_len):
        alleles = [h[i] for h in haplotypes]
        # We assume biallelic sites 'A' (ref) and 'T' (alt)
        k = alleles.count('T') # Count of alternate allele
        # Number of pairs with different alleles at this site is k * (n-k)
        pairwise_diffs_sum += k * (n_haplotypes - k)
        
    # Pi is the average number of differences per site
    # Total pairwise differences / total number of pairs
    pi = pairwise_diffs_sum / num_pairs
    return pi / seq_len


def calculate_an(n):
    """Calculates the harmonic series number a_n for Watterson's estimator."""
    if n < 2:
        return 1.0
    return sum(1.0 / i for i in range(1, n))

def run_simulation():
    """
    Simulates the described bioinformatics pipeline and calculates the bias.
    """
    # --- 1. Simulation Parameters ---
    n_haplotypes = 50  # Number of sampled haplotypes (25 individuals)
    seq_len = 10000    # Length of the genomic region
    true_n_snps = 100  # True number of Single Nucleotide Polymorphisms (SNPs)
    filter_fraction = 0.3 # A "substantial minority" (30%) of SNVs per sample are filtered

    print("--- Simulation Setup ---")
    print(f"Haplotypes: {n_haplotypes}")
    print(f"True number of SNVs (S): {true_n_snps}")
    print(f"Filter fraction per haplotype: {filter_fraction}\n")

    # --- 2. Generate 'True' Genetic Data ---
    # Start with all haplotypes being the reference ('A')
    ref_allele = 'A'
    alt_allele = 'T'
    true_haplotypes = [[ref_allele] * seq_len for _ in range(n_haplotypes)]
    
    # variants_per_haplotype[i] will store the positions of SNVs for haplotype i
    variants_per_haplotype = [[] for _ in range(n_haplotypes)]
    
    snp_positions = sorted(random.sample(range(seq_len), true_n_snps))

    for pos in snp_positions:
        # For each SNP, decide how many haplotypes carry it (must be at least 1)
        k = random.randint(1, n_haplotypes - 1)
        hap_indices_with_snp = random.sample(range(n_haplotypes), k)
        for hap_idx in hap_indices_with_snp:
            true_haplotypes[hap_idx][pos] = alt_allele
            variants_per_haplotype[hap_idx].append(pos)

    # --- 3. Calculate True Genetic Diversity Statistics ---
    s_true = true_n_snps
    an = calculate_an(n_haplotypes)
    theta_w_true = s_true / an
    pi_true = calculate_pi(true_haplotypes)

    print("--- Step 1: True Genetic Diversity (Before aaulty processing) ---")
    print(f"True Segregating Sites (S): {s_true}")
    print(f"Watterson's theta (True): {theta_w_true:.4f} = {s_true} / {an:.4f}")
    print(f"Nucleotide Diversity (pi) (True): {pi_true:.6f}\n")

    # --- 4. Simulate Filtering and Imputation ---
    # This creates the 'observed' data after the faulty pipeline
    observed_haplotypes = [[ref_allele] * seq_len for _ in range(n_haplotypes)]
    
    for hap_idx in range(n_haplotypes):
        variants = variants_per_haplotype[hap_idx]
        num_to_filter = int(len(variants) * filter_fraction)
        
        # Randomly choose which variants to 'filter out' for this haplotype
        filtered_out_variants = set(random.sample(variants, num_to_filter))
        
        # Reconstruct the observed haplotype
        # Any position not explicitly kept is imputed as the reference allele
        for pos in variants:
            if pos not in filtered_out_variants:
                observed_haplotypes[hap_idx][pos] = alt_allele

    # --- 5. Calculate Observed Genetic Diversity Statistics ---
    s_obs = 0
    for i in range(seq_len):
        alleles_at_site = [h[i] for h in observed_haplotypes]
        if len(set(alleles_at_site)) > 1:
            s_obs += 1
            
    theta_w_obs = s_obs / an
    pi_obs = calculate_pi(observed_haplotypes)

    print("--- Step 2: Observed Genetic Diversity (After faulty processing) ---")
    print("Process: For each haplotype, a fraction of its true variants were randomly deleted.")
    print("The missing sites were then imputed with the reference allele.\n")
    print(f"Observed Segregating Sites (S): {s_obs}")
    print(f"Watterson's theta (Observed): {theta_w_obs:.4f} = {s_obs} / {an:.4f}")
    print(f"Nucleotide Diversity (pi) (Observed): {pi_obs:.6f}\n")

    # --- 6. Conclusion ---
    print("--- Conclusion ---")
    print("The simulation shows that the described processing pipeline leads to:")
    if theta_w_obs < theta_w_true:
        print(f"- An underestimation of Watterson's theta ({theta_w_obs:.4f} < {theta_w_true:.4f}). It is BIASED.")
    if pi_obs < pi_true:
        print(f"- An underestimation of Nucleotide Diversity (pi) ({pi_obs:.6f} < {pi_true:.6f}). It is BIASED.")
    
    print("\nThis happens because randomly removing variant information from samples and replacing it with the reference allele artificially reduces alternate allele frequencies. This can cause some variants to be lost entirely (biasing theta) and reduces the heterozygosity at remaining sites (biasing pi).")
    print("\nTherefore, both Watterson's estimator (theta) and pi (nucleotide diversity) are biased.")

if __name__ == '__main__':
    run_simulation()