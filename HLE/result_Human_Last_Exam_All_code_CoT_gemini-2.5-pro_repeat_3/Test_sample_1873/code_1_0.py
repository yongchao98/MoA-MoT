import math

# Given parameters
gamma = 2.1  # Power-law exponent
alpha = 2.5  # Pareto distribution shape
epsilon = 0.05  # Marginal completeness tolerance
delta = 0.01  # Error probability (1 - confidence level)

# Step 1: Calculate the heterogeneity factor H
# H = (gamma - 1) / (alpha - 1)
h_numerator = gamma - 1
h_denominator = alpha - 1
H = h_numerator / h_denominator

# Step 2: Calculate the ratio of statistical error terms
stat_ratio = delta / epsilon

# Step 3: Calculate the minimum sampling ratio r
# The model assumes a direct relationship: r = H * (delta / epsilon)
r = H * stat_ratio

# Round the final result to 4 decimal places
r_rounded = round(r, 4)

# Print the final equation with all the components
print(f"Calculation of the sampling ratio r:")
print(f"Heterogeneity Factor H = (γ - 1) / (α - 1) = ({gamma} - 1) / ({alpha} - 1) = {h_numerator:.1f} / {h_denominator:.1f} = {H:.4f}")
print(f"Statistical Error Ratio = δ / ε = {delta} / {epsilon} = {stat_ratio:.1f}")
print(f"Final Sampling Ratio r = H * (δ / ε) = {H:.4f} * {stat_ratio:.1f} = {r:.4f}")
print(f"The minimum sampling ratio r required is {r_rounded}.")
