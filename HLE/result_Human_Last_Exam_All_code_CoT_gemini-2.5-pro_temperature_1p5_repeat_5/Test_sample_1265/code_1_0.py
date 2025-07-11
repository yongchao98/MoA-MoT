import math

def calculate_detection_probability(allele_frequency, sequencing_depth):
    """
    Calculates the probability of observing a rare allele at least once.

    Args:
        allele_frequency (float): The frequency of the allele in the population (e.g., 0.01 for 1%).
        sequencing_depth (int): The average number of times each base is sequenced.

    Returns:
        float: The probability of detection (0.0 to 1.0).
    """
    # Probability of NOT observing the allele in a single read
    prob_miss_single_read = 1 - allele_frequency
    
    # Probability of NOT observing the allele in ANY of the reads (at a given depth)
    prob_miss_all_reads = math.pow(prob_miss_single_read, sequencing_depth)
    
    # The probability of observing it at least once is 1 minus the probability of missing it every time
    prob_detection = 1 - prob_miss_all_reads
    
    return prob_detection

# --- Parameters ---
# Let's assume the drug resistance allele has a very low frequency of 1%
low_allele_frequency = 0.01
# Original sequencing depth
depth_40x = 40
# Proposed increased sequencing depth
depth_80x = 80

# --- Calculation ---
prob_at_40x = calculate_detection_probability(low_allele_frequency, depth_40x)
prob_at_80x = calculate_detection_probability(low_allele_frequency, depth_80x)

# --- Output ---
print(f"Goal: Detect a rare allele with a frequency of {low_allele_frequency*100}%.")
print("-" * 50)
print(f"With a sequencing depth of {depth_40x}X:")
print(f"The probability of detecting the allele at least once is {prob_at_40x:.2%}.")
print("-" * 50)
print(f"By increasing the sequencing depth to {depth_80x}X:")
print(f"The probability of detecting the allele at least once is {prob_at_80x:.2%}.")
print("-" * 50)
print("Conclusion:")
print("Increasing sequencing depth from 40X to 80X substantially improves the chance of detecting a low-frequency allele.")
print("This, combined with a thorough library preparation to ensure high data quality (minimizing false positives), makes for the most effective strategy.")
