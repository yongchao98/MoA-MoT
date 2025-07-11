import random

def simulate_allele_detection(true_frequency, depth):
    """
    Simulates sequencing reads and counts allele frequencies.

    Args:
        true_frequency (float): The actual frequency of the resistant allele (e.g., 0.02 for 2%).
        depth (int): The sequencing depth (number of reads).
    """
    resistant_allele_count = 0
    wild_type_allele_count = 0

    for _ in range(depth):
        # Each read has a random chance of being the resistant allele
        if random.random() < true_frequency:
            resistant_allele_count += 1
        else:
            wild_type_allele_count += 1

    # Calculate the observed frequency
    if depth > 0:
        observed_frequency = resistant_allele_count / depth
    else:
        observed_frequency = 0

    print(f"--- Simulation at {depth}X Depth ---")
    # Final "equation" numbers representing the counts
    print(f"Equation of Counts: {wild_type_allele_count} (Wild-Type) + {resistant_allele_count} (Resistant) = {depth} (Total)")
    print(f"Observed Frequency of Resistant Allele: {observed_frequency:.4f} (or {observed_frequency:.2%})")
    if resistant_allele_count > 0:
        print("Result: Resistant allele was detected.")
    else:
        print("Result: Resistant allele was NOT detected.")
    print("-" * 35 + "\n")


# --- Main ---
# Let's assume the drug resistance allele exists at a true frequency of 1.5% in the population.
true_resistant_allele_frequency = 0.015

print(f"Identifying a rare allele with a true population frequency of {true_resistant_allele_frequency:.1%}.\n")

# Simulate the scenario with the initial 40X sequencing depth
low_depth = 40
simulate_allele_detection(true_resistant_allele_frequency, low_depth)

# Simulate the scenario with the increased 80X sequencing depth
high_depth = 80
simulate_allele_detection(true_resistant_allele_frequency, high_depth)