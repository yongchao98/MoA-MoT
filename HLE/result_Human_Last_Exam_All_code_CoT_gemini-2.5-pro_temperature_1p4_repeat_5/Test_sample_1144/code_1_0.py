import math

# Step 1: Define the constants
# Bohr radius (r1) in meters. This is the radius of the first Bohr orbit.
bohr_radius = 5.29177e-11  # meters
# Uncertainty in position (delta_x) provided by the user, converted to meters.
delta_x = 10e-12  # meters (since 1 pm = 1e-12 m)

# Step 2: Calculate the ratio
# The ratio of momentum uncertainty (delta_p) to momentum (p) can be simplified to:
# Ratio = r1 / (2 * delta_x)
# This is derived from Heisenberg's principle (delta_p = h_bar / (2 * delta_x))
# and the momentum in the first Bohr orbit (p = h_bar / r1).
ratio = bohr_radius / (2 * delta_x)

# Step 3: Print the result and the equation
print("The ratio is calculated using the formula: Ratio = r₁ / (2 * Δx)")
print(f"Where the radius of the first Bohr orbit (r₁) = {bohr_radius:.5e} m")
print(f"And the uncertainty in position (Δx) = {delta_x:.5e} m")
print("\nFinal Equation:")
print(f"Ratio = {bohr_radius:.5f} / (2 * {delta_x:.5f})")
print(f"\nCalculated Ratio: {ratio}")
