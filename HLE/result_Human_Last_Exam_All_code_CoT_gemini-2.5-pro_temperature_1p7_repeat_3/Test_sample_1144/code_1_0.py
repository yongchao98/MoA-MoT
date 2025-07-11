import math

# Constants and given values
# r1 is the Bohr radius in meters
r1 = 5.29177e-11
# delta_x is the uncertainty in position in meters
delta_x = 10e-12

# Calculate the ratio using the simplified formula: r1 / (2 * delta_x)
ratio = r1 / (2 * delta_x)

# Print the final equation and the result
print("The ratio of the uncertainty of momentum to the momentum is given by:")
print("Ratio = (Bohr Radius) / (2 * Uncertainty in Position)")
print(f"Ratio = {r1} / (2 * {delta_x})")
print(f"The calculated ratio is: {ratio}")
