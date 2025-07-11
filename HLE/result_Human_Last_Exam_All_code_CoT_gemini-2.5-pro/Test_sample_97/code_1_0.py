import math

# Constants
K_sp = 5.3e-27
K_f = 1.1e31
K_w = 1.0e-14

# --- Step 1: Solve for [OH-] ---
# The charge balance equation leads to a quadratic equation in [OH-]^2 (let z = [OH-]^2):
# a*z^2 + b*z + c = 0
# where:
# a = 1 + K_sp * K_f
# b = -K_w
# c = -3 * K_sp

ksp_kf = K_sp * K_f
a = 1 + ksp_kf
b = -K_w
c = -3 * K_sp

# Solve the quadratic equation for z = [OH-]^2 using the quadratic formula
# z = (-b + sqrt(b^2 - 4ac)) / 2a  (we take the positive root since z must be positive)
discriminant = math.sqrt(b**2 - 4*a*c)
z = (-b + discriminant) / (2*a)

# Calculate [OH-]
oh_conc = math.sqrt(z)

# --- Step 2: Calculate concentrations of dissolved species ---
# [Al^3+] = K_sp / [OH-]^3
al_ion_conc = K_sp / (oh_conc**3)

# [Al(OH)4^-] = K_sp * K_f * [OH-]
complex_ion_conc = ksp_kf * oh_conc

# --- Step 3: Calculate total solubility (S) ---
# S = [Al^3+] + [Al(OH)4^-]
solubility = al_ion_conc + complex_ion_conc

# --- Step 4: Print the results ---
print("--- Calculation of Al(OH)3 Solubility ---")
print(f"Calculated [OH-]: {oh_conc:.3e} mol/L")
print(f"Calculated [Al^3+]: {al_ion_conc:.3e} mol/L")
print(f"Calculated [Al(OH)4^-]: {complex_ion_conc:.3e} mol/L")
print("\n--- Final Solubility Equation ---")
print(f"Solubility S = [Al^3+] + [Al(OH)4^-]")
print(f"S = {al_ion_conc:.3e} mol/L + {complex_ion_conc:.3e} mol/L")
print(f"Total Solubility (S) = {solubility:.3e} mol/L")

# Final answer in the required format
# print(f"\n<<<{solubility:.3e}>>>")