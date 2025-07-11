import math
from scipy.stats import norm

# Given parameters
confidence_level = 0.99
epsilon = 0.05  # marginal completeness tolerance
gamma = 2.1     # power-law exponent
alpha = 2.5     # Pareto distribution shape

# Step 1: Calculate the Z-score for the given confidence level
# The Z-score corresponds to the critical value for a two-tailed test.
# alpha_stat is the significance level, which is 1 - confidence_level.
alpha_stat = 1 - confidence_level
# We use ppf(1 - alpha_stat / 2) for the two-tailed Z-score.
z_score = norm.ppf(1 - alpha_stat / 2)

# Step 2: Define factors based on graph structure
# The factor for scale-free property represents the increased difficulty
# as gamma approaches 2.
gamma_factor = gamma - 2

# The factor for entity neighborhood similarity represents the increased
# difficulty for heavier-tailed similarity distributions (smaller alpha).
alpha_factor = alpha - 1

# Step 3: Calculate the minimum sampling ratio r using the derived formula
# The formula combines statistical requirements and graph complexity.
# r = (Z * ε) / ((γ - 2) * (α - 1))
r = (z_score * epsilon) / (gamma_factor * alpha_factor)

# Step 4: Print the final equation and the result
print("The calculation for the minimum sampling ratio r follows the formula:")
print("r = (Z * ε) / ((γ - 2) * (α - 1))")
print("\nWhere:")
print(f" - Confidence Level = {confidence_level} -> Z-score (Z) = {z_score:.4f}")
print(f" - Marginal Tolerance (ε) = {epsilon}")
print(f" - Power-law Exponent (γ) = {gamma}")
print(f" - Pareto Shape (α) = {alpha}")

print("\nPlugging in the values:")
print(f"r = ({z_score:.4f} * {epsilon}) / (({gamma} - 2) * ({alpha} - 1))")
print(f"r = {z_score * epsilon:.4f} / (({gamma_factor:.1f}) * ({alpha_factor:.1f}))")
print(f"r = {z_score * epsilon:.4f} / {gamma_factor * alpha_factor:.2f}")
print(f"r = {r:.4f}")

# Final Answer
print(f"\nThe minimum required sampling ratio r is: {r:.4f}")
print("<<<r=0.8586>>>")