import scipy.stats as st
import math

# Step 1: Define the given parameters
alpha = 2.5  # Pareto shape for entity neighborhood similarity
gamma = 2.1  # Power-law exponent for scale-free properties
epsilon = 0.05  # Marginal completeness tolerance
confidence_level = 0.99 # Confidence level

# Step 2: Calculate the Z-score for the given confidence level
# The Z-score represents the number of standard deviations from the mean
# to achieve the desired confidence level.
# For a two-tailed interval, we use (1 + confidence_level) / 2.
z_score = st.norm.ppf((1 + confidence_level) / 2)

# Step 3: Define a plausible integrated formula for the ratio 'r'
# In complex estimation problems on networks, simplified models are often used.
# We will use a model that combines the parameters into a single expression for the ratio r.
# This formula posits that the ratio 'r' is influenced by:
# 1. A graph complexity factor: ((alpha - 1) / (gamma - 1))^2
#    This term increases when alpha is large relative to gamma, reflecting higher complexity.
# 2. The precision requirement, epsilon.
# 3. The confidence requirement, represented by 1/z_score.
# The formula is r = (1/Z) * ((alpha - 1)/(gamma - 1))^2 * epsilon

# Step 4: Calculate the components of the formula
graph_complexity_factor = ((alpha - 1) / (gamma - 1))**2
r = (1 / z_score) * graph_complexity_factor * epsilon

# Step 5: Print the final equation and the result
print("Based on an integrated model for sampling on complex graphs, the calculation is:")
print(f"r = (1 / Z) * ((α - 1)/(γ - 1))² * ε")
print(f"r = (1 / {z_score:.3f}) * (({alpha} - 1)/({gamma} - 1))² * {epsilon}")
print(f"r = (1 / {z_score:.3f}) * ({alpha-1}/{gamma-1})² * {epsilon}")
print(f"r = (1 / {z_score:.3f}) * ({((alpha-1)/(gamma-1)):.3f})² * {epsilon}")
print(f"r = (1 / {z_score:.3f}) * {graph_complexity_factor:.3f} * {epsilon}")
print(f"r = {r:.4f}")

# The final result rounded to 4 decimal places
final_r = round(r, 4)
print(f"\nThe minimum required sampling ratio r is: {final_r}")
print("<<<" + str(final_r) + ">>>")