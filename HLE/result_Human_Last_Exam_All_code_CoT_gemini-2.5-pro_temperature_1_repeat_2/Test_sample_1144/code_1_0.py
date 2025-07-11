import math

# This script calculates the ratio of the uncertainty of an electron's momentum
# to its momentum in the first Bohr orbit.

# Define the constants and given values.
# r_1 is the radius of the first Bohr orbit (Bohr radius) in meters.
r_1 = 5.29e-11
# delta_x is the given uncertainty in position in meters (10 pm = 10e-12 m).
delta_x = 10e-12

# The ratio is derived from Heisenberg's Uncertainty Principle and the Bohr model.
# Uncertainty in momentum, delta_p = h_bar / (2 * delta_x)
# Momentum in first Bohr orbit, p = h_bar / r_1
# The ratio delta_p / p simplifies to r_1 / (2 * delta_x).
ratio = r_1 / (2 * delta_x)

# Print the final equation with the numerical values.
print(f"The equation for the ratio is r_1 / (2 * Δx)")
print(f"Plugging in the values: {r_1} / (2 * {delta_x})")

# Print the final calculated ratio.
print(f"The calculated ratio is: {ratio}")
<<<2.645>>>