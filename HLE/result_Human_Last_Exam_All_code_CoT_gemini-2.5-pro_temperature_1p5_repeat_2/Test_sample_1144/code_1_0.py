import math

# Step 1: Define the given values and constants.
# Uncertainty in position in meters (10 pm = 10e-12 m)
delta_x = 10e-12

# Radius of the first Bohr orbit in meters (Bohr radius, a0)
r1 = 5.29177e-11

# Step 2: Calculate the ratio using the simplified formula: ratio = r1 / (2 * delta_x)
# The uncertainty in momentum (delta_p) is h_bar / (2 * delta_x)
# The momentum in the first Bohr orbit (p1) is h_bar / r1
# The ratio delta_p / p1 simplifies to r1 / (2 * delta_x) as h_bar cancels out.
ratio = r1 / (2 * delta_x)

# Step 3: Print the final equation with the numbers and the result.
print("The ratio of the uncertainty of the momentum to the momentum is given by:")
print("Ratio = (Radius of first Bohr orbit) / (2 * Uncertainty of position)")
print(f"Ratio = {r1} m / (2 * {delta_x} m)")
print(f"Ratio = {ratio}")
