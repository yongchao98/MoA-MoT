#
# This script calculates the ratio of the uncertainty of an electron's momentum
# to its momentum in the first Bohr orbit.
#

# Step 1: Define the known constants.
# The Bohr radius (a0) in meters.
bohr_radius_a0 = 5.29177e-11  # meters

# The given uncertainty in position (delta_x) in meters.
# 10 pm = 10 * 10^-12 meters
uncertainty_position_dx = 10e-12  # meters

# Step 2: Calculate the ratio using the derived formula: Ratio = a0 / (2 * dx).
# The Heisenberg Uncertainty Principle gives the uncertainty in momentum (dp):
# dp = h_bar / (2 * dx)
# The momentum in the first Bohr orbit (p) is:
# p = h_bar / a0
# The ratio dp / p simplifies to a0 / (2 * dx).
ratio = bohr_radius_a0 / (2 * uncertainty_position_dx)

# Step 3: Print the final equation with the numbers used and the result.
print("The final equation for the ratio is:")
print(f"Ratio = Bohr Radius / (2 * Uncertainty in Position)")
print(f"Ratio = {bohr_radius_a0} / (2 * {uncertainty_position_dx})")
print(f"The calculated ratio is: {ratio}")