import math

# Define the physical constants and given values
# Bohr radius (a_0) in meters
a_0 = 5.29177e-11  # meters

# Given uncertainty in position (delta_x) in picometers and convert to meters
delta_x_pm = 10
delta_x_m = delta_x_pm * 1e-12 # meters

# The ratio of the uncertainty in momentum (delta_p) to the momentum (p) is derived as:
# delta_p / p = (ħ / (2 * delta_x)) / (ħ / a_0) = a_0 / (2 * delta_x)
# where ħ is the reduced Planck constant.

# Calculate the ratio
ratio = a_0 / (2 * delta_x_m)

# Print the final equation with the numbers plugged in
print("The ratio is calculated using the formula: Ratio = a_0 / (2 * Δx)")
print(f"Ratio = {a_0} / (2 * {delta_x_m})")
print("\nFinal calculated ratio:")
print(ratio)