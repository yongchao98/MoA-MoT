import math

# This script calculates the second osmotic virial coefficient from steric-only
# behavior for a typical monoclonal antibody (mAb).

# --- Step 1: Define physical constants and standard mAb parameters ---
# Avogadro's number (mol^-1)
NA = 6.022e23
# Pi
PI = math.pi
# Molar Mass (M) of a typical mAb (g/mol)
M = 150000.0
# Hydrodynamic radius (r_h) of a typical mAb in nm
r_h_nm = 5.5
# Convert radius from nanometers to centimeters for volume calculation
r_h_cm = r_h_nm * 1e-7

# --- Step 2: Calculate the molar second virial coefficient for hard spheres (B2,HS) ---
# The volume of a single mAb molecule treated as a sphere (cm^3)
v_molecule_cm3 = (4.0 / 3.0) * PI * (r_h_cm**3)
# The molar second virial coefficient (cm^3/mol, which is equivalent to mL/mol)
# B2,HS = 4 * NA * V_molecule
b2_hs_ml_per_mol = 4.0 * NA * v_molecule_cm3

# --- Step 3: Calculate the specific second virial coefficient in mL/g ---
# This quantity is B2,HS divided by the molar mass M
b22_steric_ml_per_g = b2_hs_ml_per_mol / M

# --- Step 4: Print the final result and the equation with substituted values ---
print("The second osmotic virial coefficient from steric-only behavior (B22, steric) is calculated using the hard-sphere model.")
print("Equation: B22, steric = (4 * NA * (4/3) * pi * r_h^3) / M\n")
print("Using the following values:")
print(f"  NA = {NA:.3e} mol^-1")
print(f"  pi = {PI:.5f}")
print(f"  r_h = {r_h_nm} nm = {r_h_cm:.1e} cm")
print(f"  M = {M:.1f} g/mol\n")

print("Final calculation with numbers:")
# Since r_h_cm^3 is a very small number, we show the intermediate calculation of B2,HS for clarity
print(f"B22, steric = ({b2_hs_ml_per_mol:.3e} mL/mol) / ({M:.1f} g/mol)")
# Or printing the full equation:
print(f"B22, steric = (4 * {NA:.3e} mol^-1 * (4/3) * {PI:.5f} * ({r_h_cm:.1e} cm)^3) / {M:.1f} g/mol")
print(f"\nThe calculated value is: {b22_steric_ml_per_g:.3f} mL/g")

<<<11.193>>>