import math

# --- Physical Constants ---
# Bohr radius in meters
a0 = 5.29177e-11
# Uncertainty in position in meters (10 pm = 10e-12 m)
delta_x = 10e-12

# --- Calculation ---
# The ratio of the uncertainty in momentum (delta_p) to the momentum (p) is:
# delta_p / p = (h_bar / (2 * delta_x)) / (h_bar / a0)
# This simplifies to:
# ratio = a0 / (2 * delta_x)
ratio = a0 / (2 * delta_x)

# --- Output the result ---
print("This script calculates the ratio of the uncertainty of an electron's momentum to its momentum in the first Bohr orbit.")
print("\nGiven:")
print(f"Uncertainty in position (Δx) = {delta_x:.1e} m")
print("\nConstants used:")
print(f"Bohr radius (a₀) = {a0:.5e} m")
print("\nThe formula for the ratio (Δp / p) simplifies to a₀ / (2 * Δx).")
print("\nCalculation:")
print(f"Ratio = {a0:.5f} / (2 * {delta_x})")
print(f"Ratio = {a0:.5f} / {2 * delta_x}")
print(f"Final Ratio = {ratio:.4f}")

print(f"\n<<<The final ratio is {ratio:.4f}>>>")