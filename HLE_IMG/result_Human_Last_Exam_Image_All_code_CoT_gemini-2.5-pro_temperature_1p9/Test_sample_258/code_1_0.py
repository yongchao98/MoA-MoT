import math

# Given length of the segment AB
ab_length = 18468

# The derivation shows that CP = AB * (1 + sqrt(5)) / 2
# Let's calculate the values for the equation
sqrt_5 = math.sqrt(5)
golden_ratio = (1 + sqrt_5) / 2

# Calculate the length of CP
cp_length = ab_length * golden_ratio

# Print the equation with the numerical values as requested
print(f"The calculation is based on the formula: CP = AB * (1 + \u221A5) / 2")
print(f"CP = {ab_length} * (1 + {sqrt_5}) / 2")

# Print the final result formatted to 4 decimal places
print(f"The length of CP is: {cp_length:.4f}")
