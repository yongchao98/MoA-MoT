import math

# Step 1: Define the given constants.
# Uncertainty in position (Δx) in meters. 10 pm = 10 * 10^-12 m.
delta_x = 10e-12

# Bohr radius (r₁) for the first orbit (n=1) in meters.
r1 = 5.29177e-11

# Step 2: The ratio (Δp / p) simplifies to r₁ / (2 * Δx).
# Let's calculate this ratio.
ratio = r1 / (2 * delta_x)

# Step 3: Print the final equation with the numbers and the result.
print("To find the ratio of the uncertainty of the momentum to the momentum (Δp/p), we use the formula:")
print("Ratio = r₁ / (2 * Δx)")
print("\nWhere:")
print(f"r₁ (Bohr radius) = {r1:.5e} m")
print(f"Δx (Position Uncertainty) = {delta_x:.5e} m")
print("\nPlugging in the values, the final equation is:")
# The format below shows each number in the final equation.
print(f"Ratio = {r1:.5f} / (2 * {delta_x:.5f})")
print(f"\nFinal calculated ratio:")
print(ratio)
