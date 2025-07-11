import math

# Step 1: Define the constants and given values.
# Bohr radius in meters (m)
a0 = 5.29177e-11
# Given uncertainty in position in meters (m). 10 pm = 10e-12 m
delta_x = 10e-12

# Step 2: Calculate the ratio using the simplified formula: ratio = a0 / (2 * delta_x)
# The derivation is as follows:
# Uncertainty in momentum (from Heisenberg's principle): delta_p = h_bar / (2 * delta_x)
# Momentum in the first Bohr orbit: p1 = h_bar / a0
# Ratio = delta_p / p1 = (h_bar / (2 * delta_x)) / (h_bar / a0) = a0 / (2 * delta_x)
denominator = 2 * delta_x
ratio = a0 / denominator

# Step 3: Print the final equation and the result.
print(f"The ratio is calculated as: a0 / (2 * delta_x)")
print(f"Ratio = {a0} / (2 * {delta_x})")
print(f"Ratio = {a0} / {denominator}")
print(f"The ratio of the uncertainty of the momentum to its momentum is: {ratio}")
