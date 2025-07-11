import scipy.constants

# This script calculates the ratio of the uncertainty of an electron's momentum 
# to its momentum in the first Bohr orbit.

# --- Given values and Constants ---
# The uncertainty in position is given as 10 pm.
# 1 pm = 1e-12 m
delta_x = 10e-12

# The Bohr radius (a₀) is a physical constant. We get its value from scipy.
a0 = scipy.constants.physical_constants['Bohr radius'][0]

# --- Derivation ---
print("Step 1: The uncertainty in momentum (Δp) is found using the Heisenberg Uncertainty Principle (Δx * Δp ≥ ħ/2).")
print("   For minimum uncertainty, Δp = ħ / (2 * Δx)")
print("\nStep 2: The momentum (p) of an electron in the first Bohr orbit is p = ħ / a₀, where a₀ is the Bohr radius.")
print("\nStep 3: The ratio (Δp / p) is calculated by dividing the expression from Step 1 by the expression from Step 2.")
print("   Ratio = [ħ / (2 * Δx)] / [ħ / a₀]")
print("   The reduced Planck constants (ħ) cancel, simplifying the formula to:")
print("   Ratio = a₀ / (2 * Δx)")

# --- Calculation ---
print("\n--- Calculation with values ---")
print(f"Bohr radius (a₀) = {a0} m")
print(f"Uncertainty in position (Δx) = {delta_x} m")

# Calculate the ratio
ratio = a0 / (2 * delta_x)

# --- Final Equation and Result ---
print("\nThe final equation is:")
# Output each number in the final equation as requested
print(f"{a0} / (2 * {delta_x})")

print("\nThe final calculated ratio is:")
print(ratio)
<<<2.6458860545150003>>>