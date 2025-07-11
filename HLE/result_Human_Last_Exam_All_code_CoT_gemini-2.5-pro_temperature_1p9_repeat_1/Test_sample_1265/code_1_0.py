import math
from scipy.stats import binom

def calculate_detection_probability(depth, frequency, min_reads):
    """
    Calculates the probability of seeing an allele at least `min_reads` times.
    P(k >= min_reads) = 1 - P(k < min_reads)
    """
    
    # Calculate the cumulative probability of seeing less than min_reads
    prob_lt_min = sum(binom.pmf(k, depth, frequency) for k in range(min_reads))
    
    # The probability of seeing at least min_reads is 1 minus that
    prob_ge_min = 1 - prob_lt_min
    
    return prob_ge_min, [binom.pmf(k, depth, frequency) for k in range(min_reads)]

# --- Parameters ---
allele_frequency = 0.01  # A rare allele at 1% frequency
min_required_reads = 3   # Minimum number of reads to trust a variant call
depth_1 = 40
depth_2 = 80

# --- Calculations ---
prob_1, summands_1 = calculate_detection_probability(depth_1, allele_frequency, min_required_reads)
prob_2, summands_2 = calculate_detection_probability(depth_2, allele_frequency, min_required_reads)

# --- Output Results ---
print("This script demonstrates why increasing sequencing depth is critical for finding rare alleles.")
print(f"Goal: Detect a rare allele with frequency f = {allele_frequency} by observing it at least {min_required_reads} times.\n")

print(f"--- At {depth_1}X Depth ---")
sum_str_1 = " + ".join([f"{p:.4f}" for p in summands_1])
print(f"The probability of seeing the allele on 0, 1, or 2 reads is: P(k=0) + P(k=1) + P(k=2)")
print(f"= {sum_str_1} = {sum(summands_1):.4f}")
print(f"Probability of successful detection (k >= {min_required_reads}): 1 - {sum(summands_1):.4f} = {prob_1:.4f} (or {prob_1:.2%})")

print(f"\n--- At {depth_2}X Depth ---")
sum_str_2 = " + ".join([f"{p:.4f}" for p in summands_2])
print(f"The probability of seeing the allele on 0, 1, or 2 reads is: P(k=0) + P(k=1) + P(k=2)")
print(f"= {sum_str_2} = {sum(summands_2):.4f}")
print(f"Probability of successful detection (k >= {min_required_reads}): 1 - {sum(summands_2):.4f} = {prob_2:.4f} (or {prob_2:.2%})")

print(f"\nConclusion: Doubling the sequencing depth from {depth_1}X to {depth_2}X increases the chance of detection nearly {prob_2/prob_1:.1f}-fold in this scenario.")
