import math
from scipy.stats import binom

def calculate_detection_probability(depth, allele_frequency, min_observations):
    """
    Calculates the probability of observing a rare allele at least a certain
    number of times given the sequencing depth and allele frequency.

    This uses the cumulative distribution function (CDF) of a binomial distribution.
    The probability of seeing the allele 'min_observations' or more times is:
    P(X >= k) = 1 - P(X < k) = 1 - P(X <= k-1)
    """
    # The number of trials is the sequencing depth
    n = depth
    # The probability of success on a single trial is the allele frequency
    p = allele_frequency
    # We want to find the probability of k >= min_observations
    k = min_observations

    # Calculate P(X <= k-1) using the binomial CDF
    prob_less_than_k = binom.cdf(k - 1, n, p)

    # The desired probability is 1 minus the above
    prob_at_least_k = 1 - prob_less_than_k
    
    return prob_at_least_k

# --- Parameters of our scenario ---
# Let's assume the drug resistance allele is present at a low frequency
allele_freq = 0.01  # 1% frequency

# We need to see the allele more than once to be confident it's not a sequencing error
min_reads_for_detection = 2 

# Original sequencing depth
depth_1 = 40

# Increased sequencing depth
depth_2 = 80

# --- Calculations ---
prob_at_40X = calculate_detection_probability(depth_1, allele_freq, min_reads_for_detection)
prob_at_80X = calculate_detection_probability(depth_2, allele_freq, min_reads_for_detection)

# --- Output the results ---
print("This script demonstrates the benefit of increasing sequencing depth for detecting rare alleles.\n")
print(f"Scenario: Identifying a rare allele with a frequency of {allele_freq*100}%.")
print(f"Condition for detection: Observing the allele on at least {min_reads_for_detection} separate reads.\n")

print("--- Analysis for 40X Depth ---")
print(f"The equation calculates the probability of detection at a depth of {depth_1}X.")
print(f"P(detection) = 1 - P(observations < {min_reads_for_detection}) for a binomial distribution with n={depth_1} and p={allele_freq}")
print(f"Result: The probability of detecting the allele is {prob_at_40X:.2%}\n")

print("--- Analysis for 80X Depth ---")
print(f"The equation calculates the probability of detection at a depth of {depth_2}X.")
print(f"P(detection) = 1 - P(observations < {min_reads_for_detection}) for a binomial distribution with n={depth_2} and p={allele_freq}")
print(f"Result: The probability of detecting the allele is {prob_at_80X:.2%}\n")

print(f"Conclusion: Increasing the depth from {depth_1}X to {depth_2}X increases the probability of confident detection from {prob_at_40X:.2%} to {prob_at_80X:.2%}, making detection more than {prob_at_80X/prob_at_40X:.1f} times as likely.")