import math

def calculate_detection_probability(depth, freq, min_hits):
    """
    Calculates the probability of observing an allele at least 'min_hits' times,
    given a certain sequencing 'depth' and allele 'freq'.
    This uses the cumulative binomial distribution formula:
    P(X >= k) = 1 - P(X < k)
    """
    prob_of_seeing_fewer_than_min_hits = 0.0
    # Calculate the sum of probabilities for observing the allele 0, 1, ..., min_hits-1 times
    for i in range(min_hits):
        try:
            # Binomial probability: C(n, k) * p^k * (1-p)^(n-k)
            prob_i = math.comb(depth, i) * (freq**i) * ((1 - freq)**(depth - i))
            prob_of_seeing_fewer_than_min_hits += prob_i
        except (ValueError, OverflowError):
            # This can happen for very small probabilities; it's effectively 0.
            pass
    
    # The probability of successful detection is 1 minus the probability of missing it
    prob_detection = 1 - prob_of_seeing_fewer_than_min_hits
    return prob_detection

# --- Parameters based on the problem ---
allele_frequency = 0.01  # A very low frequency allele at 1%
initial_depth = 40
increased_depth = 80
# We need to see the allele at least 3 times to be confident it's not a sequencing error
min_reads_for_confidence = 3 

# --- Perform Calculations ---
prob_at_40x = calculate_detection_probability(initial_depth, allele_frequency, min_reads_for_confidence)
prob_at_80x = calculate_detection_probability(increased_depth, allele_frequency, min_reads_for_confidence)

# --- Print the results and the "equation" numbers ---
print("This script demonstrates why increasing sequencing depth is crucial for detecting rare alleles.")
print(f"Goal: Confidently detect an allele with a frequency of {allele_frequency*100}%.")
print(f"Confidence threshold: Must observe the allele on at least {min_reads_for_confidence} reads.\n")

print("--- Case 1: Initial Depth (40X) ---")
print(f"The equation uses these numbers: Depth (n) = {initial_depth}, Allele Frequency (p) = {allele_frequency}, Minimum Hits (k) = {min_reads_for_confidence}")
print(f"Probability of confident detection at {initial_depth}X depth: {prob_at_40x:.4f} (or {prob_at_40x*100:.2f}%)")
print("-" * 20)

print("\n--- Case 2: Increased Depth (80X) ---")
print(f"The equation uses these numbers: Depth (n) = {increased_depth}, Allele Frequency (p) = {allele_frequency}, Minimum Hits (k) = {min_reads_for_confidence}")
print(f"Probability of confident detection at {increased_depth}X depth: {prob_at_80x:.4f} (or {prob_at_80x*100:.2f}%)")
print("-" * 20)

print(f"\nConclusion: By increasing sequencing depth from {initial_depth}X to {increased_depth}X, the chance of confidently detecting the rare allele improves dramatically.")