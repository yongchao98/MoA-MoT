import math

def calculate_detection_probability(depth, allele_frequency, min_reads):
    """
    Calculates the probability of observing at least 'min_reads' for a rare allele.
    This uses the cumulative distribution function of the binomial distribution.
    P(X >= k) = 1 - P(X < k) = 1 - sum_{i=0}^{k-1} P(X=i)
    """
    prob_less_than_min = 0.0
    for i in range(min_reads):
        # Binomial probability: P(X=i) = C(n, k) * p^k * (1-p)^(n-k)
        try:
            prob_i = (math.comb(depth, i) *
                      (allele_frequency ** i) *
                      ((1 - allele_frequency) ** (depth - i)))
            prob_less_than_min += prob_i
        except (ValueError, OverflowError):
            # Handle cases where numbers are too large for standard floats
            return float('nan')
            
    return 1 - prob_less_than_min

# --- Parameters ---
# Allele frequency for a "very low frequency" variant
f = 0.01  # 1%
# Minimum number of reads required to confidently call a variant
k = 3
# Sequencing depths from the problem
depth1 = 40
depth2 = 80

# --- Calculation ---
prob_at_40x = calculate_detection_probability(depth1, f, k)
prob_at_80x = calculate_detection_probability(depth2, f, k)

# --- Output ---
print("This script calculates the probability of detecting a rare allele.")
print(f"Assuming an allele frequency of {f*100}% and requiring at least {k} reads for detection:\n")

# Output for 40X depth
print(f"At {depth1}X depth:")
# The prompt asks to output each number in the final equation.
# The equation is P(X >= k) = 1 - sum(P(X=i) for i=0..k-1)
print(f"The probability of detection is calculated as: 1 - sum(C({depth1}, i) * {f}^i * (1-{f})^({depth1}-i) for i=0 to {k-1})")
print(f"Result: {prob_at_40x:.4f}\n")

# Output for 80X depth
print(f"At {depth2}X depth:")
print(f"The probability of detection is calculated as: 1 - sum(C({depth2}, i) * {f}^i * (1-{f})^({depth2}-i) for i=0 to {k-1})")
print(f"Result: {prob_at_80x:.4f}\n")

print(f"By increasing depth from {depth1}X to {depth2}X, the probability of detecting the rare allele increases significantly.")
