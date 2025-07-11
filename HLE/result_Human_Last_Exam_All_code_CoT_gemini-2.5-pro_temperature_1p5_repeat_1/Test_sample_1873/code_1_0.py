import scipy.stats
import math

# Step 1: Define the given parameters
gamma = 2.1  # Power-law exponent for scale-free property
alpha = 2.5  # Shape parameter for truncated Pareto distribution of neighborhood similarity
epsilon = 0.05  # Marginal completeness tolerance
confidence_level = 0.99  # Confidence level

# Calculate the Z-score for the given confidence level
# The significance level (delta) is 1 - confidence_level
delta = 1 - confidence_level
# For a two-tailed interval, we use delta / 2
z_score = scipy.stats.norm.ppf(1 - delta / 2)

# Step 2: Formulate and calculate the final ratio r
# The model for r combines a design effect factor and a statistical requirement factor.
# Design Effect Factor (Deff_factor) based on gamma and alpha
# This factor models the increased variance due to the graph's complex structure.
deff_factor_gamma = (gamma - 1) / (gamma - 2)
deff_factor_alpha = (alpha - 1) / (alpha - 2)
deff_factor = deff_factor_gamma + deff_factor_alpha

# Statistical Requirement Factor (Stat_factor) based on epsilon and Z-score
stat_factor = (epsilon / z_score)**2

# The final formula for the ratio r
r = deff_factor * stat_factor

# Step 3: Print the intermediate values and the final result
print(f"Given parameters:")
print(f"Power-law exponent (γ): {gamma}")
print(f"Pareto shape (α): {alpha}")
print(f"Tolerance (ε): {epsilon}")
print(f"Confidence Level: {confidence_level}\n")

print(f"Derived values:")
print(f"Z-score for {confidence_level*100}% confidence: {z_score:.4f}")
print(f"Design effect contribution from γ: (γ-1)/(γ-2) = {deff_factor_gamma:.4f}")
print(f"Design effect contribution from α: (α-1)/(α-2) = {deff_factor_alpha:.4f}")
print(f"Combined design effect factor: {deff_factor:.4f}\n")

print(f"Final Equation for r:")
print(f"r = ( (γ-1)/(γ-2) + (α-1)/(α-2) ) * (ε/Z)²")
print(f"r = ( ({gamma}-1)/({gamma}-2) + ({alpha}-1)/({alpha}-2) ) * ({epsilon}/{z_score:.4f})²")
print(f"r = ( {deff_factor_gamma:.4f} + {deff_factor_alpha:.4f} ) * ({epsilon/z_score:.4f})²")
print(f"r = {deff_factor:.4f} * {stat_factor:.8f}")
print(f"r = {r:.8f}\n")

# Output the final answer rounded to 4 decimal places
final_r = round(r, 4)
print(f"The minimum required sampling ratio r is: {final_r}")
<<<0.0053>>>