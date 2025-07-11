import random
import math

def run_simulation():
    """
    Simulates the effect of per-sample filtering and reference imputation
    on Watterson's theta and nucleotide diversity (pi).
    """
    # --- 1. Parameters ---
    n_haplotypes = 100  # Number of sequences (samples)
    seq_len = 50000   # Length of the sequence region
    true_n_seg_sites = 200 # True number of polymorphic sites
    p_filter = 0.25   # 25% chance a true variant call is low-quality and filtered

    # --- 2. Generate "True" Data ---
    # Start with all reference alleles ('A')
    true_haplotypes = [['A'] * seq_len for _ in range(n_haplotypes)]

    # Select sites to be polymorphic
    poly_sites_indices = random.sample(range(seq_len), true_n_seg_sites)

    # Introduce variants ('T') at these sites
    for site_idx in poly_sites_indices:
        # For each polymorphic site, decide how many haplotypes carry the variant
        # Let's make 1 to 50% of haplotypes carry the variant
        num_variants = random.randint(1, n_haplotypes // 2)
        variant_haplotype_indices = random.sample(range(n_haplotypes), num_variants)
        for hap_idx in variant_haplotype_indices:
            true_haplotypes[hap_idx][site_idx] = 'T'

    # --- 3. Simulate Filtering and Imputation to create "Processed" Data ---
    processed_haplotypes = [row[:] for row in true_haplotypes] # Deep copy
    for i in range(n_haplotypes):
        for j in range(seq_len):
            # If a true variant is found
            if true_haplotypes[i][j] == 'T':
                # Simulate a random chance of it being filtered
                if random.random() < p_filter:
                    # Impute with the reference allele
                    processed_haplotypes[i][j] = 'A'

    # --- 4. Define Calculation Functions ---
    def calculate_s(haplotypes):
        """Calculates the number of segregating sites (S)."""
        s_count = 0
        num_sites = len(haplotypes[0])
        for j in range(num_sites):
            alleles = set(haplotypes[i][j] for i in range(len(haplotypes)))
            if len(alleles) > 1:
                s_count += 1
        return s_count

    def calculate_pi_sum(haplotypes):
        """Calculates the sum of pairwise nucleotide diversity across all sites."""
        n = len(haplotypes)
        num_sites = len(haplotypes[0])
        total_pi_sum = 0
        for j in range(num_sites):
            # Count variant allele frequency 'p'
            variant_count = sum(1 for i in range(n) if haplotypes[i][j] == 'T')
            p = variant_count / n
            # Site's contribution to pi is 2*p*(1-p)
            total_pi_sum += 2 * p * (1 - p)
        return total_pi_sum

    # --- 5. Perform Calculations and Print Results ---
    # a_n is the (n-1)th harmonic number
    a_n = sum(1.0/i for i in range(1, n_haplotypes))

    print("--- Analysis of 'True' Data (Before Filtering) ---")
    s_true = calculate_s(true_haplotypes)
    pi_sum_true = calculate_pi_sum(true_haplotypes)
    theta_w_true = s_true / a_n
    pi_true = pi_sum_true / seq_len
    
    print(f"Number of Segregating Sites (S): {s_true}")
    print(f"Harmonic number (a_n for n={n_haplotypes}): {a_n:.4f}")
    print(f"Watterson's theta (S / a_n): {s_true} / {a_n:.4f} = {theta_w_true:.4f}\n")
    
    print(f"Sum of pairwise diversity (Pi_sum): {pi_sum_true:.4f}")
    print(f"Sequence Length (L): {seq_len}")
    print(f"Nucleotide Diversity (pi = Pi_sum / L): {pi_sum_true:.4f} / {seq_len} = {pi_true:.6f}\n")

    print("--- Analysis of 'Processed' Data (After Filtering) ---")
    s_processed = calculate_s(processed_haplotypes)
    pi_sum_processed = calculate_pi_sum(processed_haplotypes)
    theta_w_processed = s_processed / a_n
    pi_processed = pi_sum_processed / seq_len

    print(f"Number of Segregating Sites (S): {s_processed}")
    print(f"Harmonic number (a_n for n={n_haplotypes}): {a_n:.4f}")
    print(f"Watterson's theta (S / a_n): {s_processed} / {a_n:.4f} = {theta_w_processed:.4f}\n")

    print(f"Sum of pairwise diversity (Pi_sum): {pi_sum_processed:.4f}")
    print(f"Sequence Length (L): {seq_len}")
    print(f"Nucleotide Diversity (pi = Pi_sum / L): {pi_sum_processed:.4f} / {seq_len} = {pi_processed:.6f}\n")
    
    print("--- Conclusion ---")
    print(f"The number of segregating sites is underestimated ({s_processed} < {s_true}), biasing Watterson's theta downwards.")
    print(f"The sum of pairwise diversity is also underestimated ({pi_sum_processed:.4f} < {pi_sum_true:.4f}), biasing nucleotide diversity (pi) downwards.")
    print("Therefore, both Watterson's estimator (theta) and pi (nucleotide diversity) are biased.")

if __name__ == '__main__':
    run_simulation()