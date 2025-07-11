import math

def calculate_detection_probability(allele_frequency, depth):
    """
    Calculates the probability of failing to detect a rare allele.
    This is the probability of observing zero reads for that allele,
    which follows a binomial distribution P(X=0) = (1-p)^n.
    """
    prob_miss = (1 - allele_frequency) ** depth
    return prob_miss

# --- Parameters ---
# We are looking for a very low frequency allele, e.g., at 1% frequency.
allele_freq = 0.01

# We compare the two sequencing depths mentioned in the problem.
depth1 = 40
depth2 = 80

# --- Calculations ---
# Calculate the probability of missing the allele at 40X depth
prob_miss_d1 = calculate_detection_probability(allele_freq, depth1)

# Calculate the probability of missing the allele at 80X depth
prob_miss_d2 = calculate_detection_probability(allele_freq, depth2)

# --- Output ---
print(f"Goal: To detect a rare allele with a frequency of {allele_freq*100}%.")
print("-" * 60)
print(f"Analysis for Sequencing Depth: {depth1}X")
print("The probability of NOT detecting the allele (i.e., getting 0 reads)")
print("is calculated using the binomial probability formula:")
print(f"P(miss) = (1 - allele_frequency)^depth")
print(f"P(miss) = (1 - {allele_freq})^{depth1} = {prob_miss_d1:.4f}")
print(f"This means there is a {prob_miss_d1*100:.2f}% chance of completely missing the allele.")
print("-" * 60)
print(f"Analysis for Sequencing Depth: {depth2}X")
print(f"By increasing the depth to {depth2}X, the calculation is:")
print(f"P(miss) = (1 - {allele_freq})^{depth2} = {prob_miss_d2:.4f}")
print(f"The chance of missing the allele drops to {prob_miss_d2*100:.2f}%.")
print("-" * 60)
print("Conclusion:")
print(f"Increasing sequencing depth from {depth1}X to {depth2}X is crucial.")
print("It significantly improves the chance of sampling the rare allele,")
print("which is a necessary first step for its identification.")
