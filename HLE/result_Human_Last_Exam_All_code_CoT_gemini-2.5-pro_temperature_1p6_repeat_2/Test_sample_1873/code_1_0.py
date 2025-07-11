import scipy.stats
import math

# Step 1: Define the given parameters
epsilon = 0.05  # Marginal completeness tolerance
confidence_level = 0.99
gamma = 2.1  # Power-law exponent
alpha = 2.5  # Pareto distribution shape

# Step 2: Calculate the Z-score for the given confidence level
# For a two-sided interval, alpha_sig is the significance level
alpha_sig = 1 - confidence_level
# The Z-score corresponds to the upper tail probability of alpha_sig / 2
z_score = scipy.stats.norm.ppf(1 - alpha_sig / 2)

# Step 3: Define the terms for the formula
# The core statistical requirement term is (Z/ε)^2
statistical_term_sq = (z_score / epsilon) ** 2

# The complexity factor from graph properties.
# Smaller gamma and alpha indicate higher heterogeneity, leading to a smaller
# complexity factor, which in turn increases the required sampling ratio.
complexity_factor = (gamma - 1) * (alpha - 1)

# Step 4: Apply the formula to calculate the sampling ratio r
# The formula is r = 1 / (statistical_term * complexity_factor)
# This model assumes the required ratio 'r' is inversely proportional
# to an effective sample size adjusted for graph complexity.
denominator = statistical_term_sq * complexity_factor
r = 1 / denominator

# Step 5: Print the results step-by-step and the final answer
print(f"Given parameters:")
print(f"Marginal tolerance (ε): {epsilon}")
print(f"Confidence level: {confidence_level}")
print(f"Power-law exponent (γ): {gamma}")
print(f"Pareto shape (α): {alpha}")
print("-" * 30)

print(f"Intermediate Calculations:")
print(f"Z-score for {confidence_level*100}% confidence: {z_score:.4f}")
print(f"Statistical requirement term (Z/ε)²: {statistical_term_sq:.4f}")
print(f"Graph complexity factor (γ-1)*(α-1): {complexity_factor:.4f}")
print("-" * 30)

print("Final Calculation of the sampling ratio r:")
# We display the full equation with the calculated numbers
print(f"r = 1 / ( ( {z_score:.4f} / {epsilon} )**2 * ( ( {gamma} - 1 ) * ( {alpha} - 1 ) ) )")
print(f"r = 1 / ( {statistical_term_sq:.4f} * {complexity_factor:.4f} )")
print(f"r = 1 / ( {denominator:.4f} )")
print(f"Calculated ratio r: {r:.10f}")
print("-" * 30)

# Final answer rounded to 4 decimal places
final_r = round(r, 4)
print(f"The minimum required sampling ratio r, rounded to 4 decimal places, is: {final_r}")

# Return final answer in the specified format
# print(f"<<<{final_r}>>>")