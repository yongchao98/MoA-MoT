import math

def analyze_detection_probability(allele_frequency, sequencing_depth):
    """
    Calculates the expected number of reads for a rare allele and the probability
    of observing it at least once.
    """
    # The expected number of reads is the depth multiplied by the frequency.
    expected_reads = sequencing_depth * allele_frequency

    # The probability of NOT seeing the allele in a single read is (1 - frequency).
    # The probability of NOT seeing it in ANY of the reads is (1 - frequency)^depth.
    prob_of_not_seeing = (1 - allele_frequency) ** sequencing_depth

    # The probability of seeing it at least once is 1 minus the probability of not seeing it.
    prob_of_seeing_at_least_once = 1 - prob_of_not_seeing

    return expected_reads, prob_of_seeing_at_least_once

# --- Parameters ---
# Let's assume a drug resistance allele is at a "very low" frequency of 0.5%.
rare_allele_frequency = 0.005

# The two sequencing depths we are comparing.
low_depth = 40
high_depth = 80

# --- Calculations ---
expected_reads_low, prob_low = analyze_detection_probability(rare_allele_frequency, low_depth)
expected_reads_high, prob_high = analyze_detection_probability(rare_allele_frequency, high_depth)

# --- Output ---
print(f"Analysis for detecting a rare allele at {rare_allele_frequency * 100}% frequency.")
print("-" * 60)

print(f"\nCase 1: Low Sequencing Depth ({low_depth}X)")
print(f"The expected number of reads supporting the rare allele is calculated as:")
print(f"Equation: {low_depth} (depth) * {rare_allele_frequency} (frequency) = {expected_reads_low:.2f} (expected reads)")
print(f"With only {expected_reads_low:.2f} expected reads, the allele is likely missed or indistinguishable from noise.")
print(f"The statistical probability of observing this allele at least once is only {prob_low:.2%}.")


print(f"\nCase 2: Increased Sequencing Depth ({high_depth}X)")
print(f"By increasing the depth, the expected number of reads becomes:")
print(f"Equation: {high_depth} (depth) * {rare_allele_frequency} (frequency) = {expected_reads_high:.2f} (expected reads)")
print(f"Doubling the depth doubles the expected reads to {expected_reads_high:.2f}.")
print(f"The statistical probability of observing this allele at least once increases to {prob_high:.2%}.")

print("-" * 60)
print("\nConclusion: Increasing the sequencing depth from 40X to 80X significantly improves the theoretical chance of detecting the rare allele, confirming it's a critical step.")