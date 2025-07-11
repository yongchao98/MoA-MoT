import random

def simulate_sequencing(depth, rare_allele_freq, error_rate):
    """
    Simulates sequencing a single genomic position to a given depth.

    Args:
        depth (int): The sequencing depth (e.g., 40 for 40X).
        rare_allele_freq (float): The frequency of the rare allele in the population.
        error_rate (float): The per-base sequencing error rate.

    Returns:
        dict: A dictionary with counts of each observed base ('A', 'T', 'C', 'G').
    """
    ref_allele = 'A'
    rare_allele = 'T'
    other_bases = ['C', 'G']

    counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0}

    for _ in range(depth):
        # Step 1: Determine the true base from the population's allele frequency
        if random.random() < rare_allele_freq:
            true_base = rare_allele
        else:
            true_base = ref_allele

        # Step 2: Simulate a potential sequencing error
        if random.random() < error_rate:
            # An error occurred, pick a random base that is not the true base
            read_base = random.choice([b for b in ['A', 'T', 'C', 'G'] if b != true_base])
        else:
            # No error occurred, the read base is the true base
            read_base = true_base

        counts[read_base] += 1

    return counts

# --- Simulation Parameters ---
# Set a seed for reproducibility
random.seed(42)

# Depths from the problem description
depth_1 = 40
depth_2 = 80

# A hypothetical low frequency for the rare allele 'T'
rare_allele_frequency = 0.025 # 2.5%

# A typical sequencing error rate for a high-accuracy platform
sequencing_error_rate = 0.01 # 1%

# --- Run Simulations ---
counts_40x = simulate_sequencing(depth_1, rare_allele_frequency, sequencing_error_rate)
counts_80x = simulate_sequencing(depth_2, rare_allele_frequency, sequencing_error_rate)

# --- Print Results ---
print(f"Goal: Detect a rare allele 'T' with a true frequency of {rare_allele_frequency*100}%")
print(f"Background sequencing error rate: {sequencing_error_rate*100}%")
print("-" * 65)

# Scenario 1: 40X Depth
print(f"Scenario 1: Sequencing at {depth_1}X depth")
num_A_40x = counts_40x.get('A', 0)
num_T_40x = counts_40x.get('T', 0)
num_err_40x = counts_40x.get('C', 0) + counts_40x.get('G', 0)
print(f"Expected 'T' reads = {depth_1} * {rare_allele_frequency} = {depth_1 * rare_allele_frequency}")
print(f"Observed reads -> A: {num_A_40x}, T: {num_T_40x}, Errors (C/G): {num_err_40x}")
print("Observation: At 40X, the rare allele 'T' is observed only once. This count is very low and could be mistaken for a random sequencing error.")

print("\n" + "-" * 65)

# Scenario 2: 80X Depth
print(f"Scenario 2: Increasing depth to {depth_2}X")
num_A_80x = counts_80x.get('A', 0)
num_T_80x = counts_80x.get('T', 0)
num_err_80x = counts_80x.get('C', 0) + counts_80x.get('G', 0)
print(f"Expected 'T' reads = {depth_2} * {rare_allele_frequency} = {depth_2 * rare_allele_frequency}")
print(f"Observed reads -> A: {num_A_80x}, T: {num_T_80x}, Errors (C/G): {num_err_80x}")
print("Observation: At 80X, the rare allele 'T' is observed 3 times. This provides stronger evidence that 'T' is a true low-frequency variant, not just noise.")
