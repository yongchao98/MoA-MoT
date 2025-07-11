import numpy as np

def run_simulation():
    """
    Simulates the bioinformatics pipeline and demonstrates its biasing effect
    on theta and pi calculations.
    """
    # --- Simulation Parameters ---
    N_HAPLOTYPES = 50  # Number of sequences (e.g., from 25 diploid samples)
    SEQ_LENGTH = 20000
    TRUE_N_SITES = 100 # Number of true segregating sites
    FILTER_PROB = 0.3  # A "substantial minority" (30%) of variants are filtered per sample

    print("--- Plan ---")
    print("1. Create a 'true' set of sequences with known diversity.")
    print("2. Calculate the 'true' theta and pi from this data.")
    print("3. Simulate the pipeline: filter variants and impute with reference.")
    print("4. Calculate the 'observed' theta and pi from the processed data.")
    print("5. Compare the 'true' and 'observed' values to show the bias.\n")
    
    # --- Step 1: Create TRUE dataset ---
    REF_ALLELE = 'A'
    ALT_ALLELE = 'T'
    true_haplotypes = np.full((N_HAPLOTYPES, SEQ_LENGTH), REF_ALLELE, dtype='<U1')

    # Randomly select polymorphic sites
    polymorphic_indices = np.random.choice(range(SEQ_LENGTH), TRUE_N_SITES, replace=False)

    # At each polymorphic site, introduce the ALT allele in a random number of haplotypes
    for site_idx in polymorphic_indices:
        n_alt_alleles = np.random.randint(1, N_HAPLOTYPES // 2)  # Ensure it's a minor allele
        alt_haplotype_indices = np.random.choice(range(N_HAPLOTYPES), n_alt_alleles, replace=False)
        true_haplotypes[alt_haplotype_indices, site_idx] = ALT_ALLELE

    # --- Step 2: Calculate TRUE theta and pi ---
    print("--- Calculating 'True' Diversity Metrics ---")

    # Watterson's Theta
    S_true = np.sum([len(np.unique(true_haplotypes[:, i])) > 1 for i in range(SEQ_LENGTH)])
    a_n = sum(1.0 / i for i in range(1, N_HAPLOTYPES))
    theta_w_true = S_true / a_n
    
    # Nucleotide Diversity (pi)
    num_pairs = N_HAPLOTYPES * (N_HAPLOTYPES - 1) / 2
    total_diffs_true = 0
    for i in range(N_HAPLOTYPES):
        for j in range(i + 1, N_HAPLOTYPES):
            total_diffs_true += np.sum(true_haplotypes[i] != true_haplotypes[j])
    pi_true = total_diffs_true / num_pairs

    print(f"True Number of Segregating Sites (S_true): {S_true}")
    print(f"Equation for theta_W: S_true / a_n")
    print(f"Calculation: {S_true} / {a_n:.4f} = {theta_w_true:.4f} (Theta for the whole sequence)")
    print("\nEquation for pi: (Sum of pairwise differences) / (Number of pairs)")
    print(f"Calculation: {total_diffs_true} / {int(num_pairs)} = {pi_true:.4f} (Pi, avg diffs per sequence pair)")


    # --- Step 3: Simulate Filtering and Imputation ---
    observed_haplotypes = np.copy(true_haplotypes)
    for i in range(N_HAPLOTYPES):
        for j in range(SEQ_LENGTH):
            # If this is a variant site for this haplotype...
            if observed_haplotypes[i, j] != REF_ALLELE:
                # ...it has a chance of being filtered and imputed with the reference
                if np.random.rand() < FILTER_PROB:
                    observed_haplotypes[i, j] = REF_ALLELE

    # --- Step 4: Calculate OBSERVED theta and pi ---
    print("\n--- Calculating 'Observed' Diversity Metrics (After Pipeline) ---")
    
    # Watterson's Theta
    S_obs = np.sum([len(np.unique(observed_haplotypes[:, i])) > 1 for i in range(SEQ_LENGTH)])
    theta_w_obs = S_obs / a_n
    
    # Nucleotide Diversity (pi)
    total_diffs_obs = 0
    for i in range(N_HAPLOTYPES):
        for j in range(i + 1, N_HAPLOTYPES):
            total_diffs_obs += np.sum(observed_haplotypes[i] != observed_haplotypes[j])
    pi_obs = total_diffs_obs / num_pairs

    print(f"Observed Number of Segregating Sites (S_obs): {S_obs}")
    print(f"Equation for theta_W: S_obs / a_n")
    print(f"Calculation: {S_obs} / {a_n:.4f} = {theta_w_obs:.4f} (Theta for the whole sequence)")
    print("\nEquation for pi: (Sum of pairwise differences) / (Number of pairs)")
    print(f"Calculation: {total_diffs_obs} / {int(num_pairs)} = {pi_obs:.4f} (Pi, avg diffs per sequence pair)")
    
    # --- Step 5: Conclusion ---
    print("\n--- Final Comparison and Conclusion ---")
    print(f"Theta_W | True: {theta_w_true:.4f} | Observed: {theta_w_obs:.4f} -> Result is BIASED")
    print(f"Pi      | True: {pi_true:.4f} | Observed: {pi_obs:.4f} -> Result is BIASED")
    print("\nBecause the observed values for both estimators are systematically lower than the true values,")
    print("both Watterson's estimator (theta) and pi (nucleotide diversity) are biased.")

run_simulation()