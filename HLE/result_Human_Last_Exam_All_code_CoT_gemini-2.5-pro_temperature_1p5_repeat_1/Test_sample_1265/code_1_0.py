import math

def binomial_pmf(k, n, p):
    """
    Calculates the Probability Mass Function (PMF) for a binomial distribution.
    This gives the probability of getting exactly 'k' successes in 'n' trials.
    P(X=k) = C(n, k) * p^k * (1-p)^(n-k)
    """
    if k < 0 or k > n:
        return 0
    # math.comb(n, k) calculates "n choose k"
    return math.comb(n, k) * (p ** k) * ((1 - p) ** (n - k))

def calculate_and_print_probabilities(depth, frequency, min_reads):
    """
    Calculates and prints the probability of detecting a low-frequency allele.
    """
    print(f"--- Analysis for {depth}X Depth ---")
    print(f"Parameters: N (depth) = {depth}, p (frequency) = {frequency}, min_reads_required = {min_reads}")

    # Calculate probability of observing 0, 1, or 2 reads (i.e., not detecting)
    prob_not_detecting = 0.0
    print("Equation for not detecting (seeing fewer than required reads):")
    print(f"P(not detect) = P(X=0) + P(X=1) + ... + P(X={min_reads-1})")

    # Show calculation for each term
    for k in range(min_reads):
        prob_k = binomial_pmf(k, depth, frequency)
        prob_not_detecting += prob_k
        print(f"P(X={k}) = C({depth}, {k}) * {frequency}^{k} * (1-{frequency})^({depth}-{k}) = {prob_k:.4f}")

    # The probability of detection is 1 minus the probability of not detecting
    prob_detecting = 1 - prob_not_detecting

    print(f"\nTotal P(not detect) = {prob_not_detecting:.4f}")
    print(f"Final P(detect) = 1 - {prob_not_detecting:.4f} = {prob_detecting:.4f}")
    print(f"Result: The probability of detection at {depth}X is {prob_detecting*100:.2f}%\n")

# --- Parameters based on the problem ---
# A very low frequency for a drug resistance allele
allele_frequency = 0.01  # 1%
# Minimum number of reads with the allele for a confident call
min_reads_required = 3
# Sequencing depths from the problem description
depth_initial = 40
depth_increased = 80

# --- Main Execution ---
print("This script demonstrates why increasing sequencing depth is critical for finding rare alleles.\n")
calculate_and_print_probabilities(depth_initial, allele_frequency, min_reads_required)
calculate_and_print_probabilities(depth_increased, allele_frequency, min_reads_required)

print("Conclusion: As shown, doubling sequencing depth from 40X to 80X significantly increases the probability of detecting a rare allele.")
print("This, combined with superior lab methods and analysis tools as described in option E, provides the most reliable strategy.")

<<<E>>>