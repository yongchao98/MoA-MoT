import math

# This script calculates the ratio of the uncertainty of an electron's momentum
# to its momentum in the first Bohr orbit.

# 1. Define the given values and constants.
# The uncertainty in position is given as 10 picometers (pm).
# 1 pm = 1e-12 meters.
delta_x = 10e-12  # in meters

# The radius of the first Bohr orbit (Bohr radius, a₀) is a physical constant.
r1 = 5.29177e-11 # in meters

# 2. Calculate the ratio using the simplified formula derived from physical principles.
# The ratio (Δp / p) simplifies to r₁ / (2 * Δx).
ratio = r1 / (2 * delta_x)

# 3. Print the explanation, the final equation with values, and the result.
print("The ratio of the uncertainty of momentum (Δp) to the momentum (p) is calculated using the formula:")
print("Ratio = r₁ / (2 * Δx)\n")

print("Given values:")
print(f"Bohr radius (r₁): {r1} m")
print(f"Uncertainty in position (Δx): {delta_x} m\n")

print("Substituting the values into the equation:")
print(f"Ratio = {r1} / (2 * {delta_x})")
print(f"Ratio = {ratio:.4f}")

print(f"\n<<<{ratio:.4f}>>>")