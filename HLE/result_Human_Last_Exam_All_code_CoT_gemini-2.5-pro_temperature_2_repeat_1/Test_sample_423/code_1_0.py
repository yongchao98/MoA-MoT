# Part 1: Analysis of B
# Given mass fractions
mass_C_percent = 50.0
mass_H_percent = 10.0
mass_O_percent = 40.0

# Molar masses
M_C = 12.011
M_H = 1.008
M_O = 15.999

# Step 1: Calculate moles assuming 100g of substance B
moles_C_B = mass_C_percent / M_C
moles_H_B = mass_H_percent / M_H
moles_O_B = mass_O_percent / M_O

# Step 2: Find the simplest ratio by dividing by the smallest number of moles
smallest_moles = min(moles_C_B, moles_H_B, moles_O_B)
ratio_C = moles_C_B / smallest_moles
ratio_H = moles_H_B / smallest_moles
ratio_O = moles_O_B / smallest_moles

# Step 3: Convert to nearest whole number ratio
# Multiply by 3 to clear the fraction for Carbon (1.66 * 3 = 5)
emp_C = round(ratio_C * 3)
emp_H = round(ratio_H * 3)
emp_O = round(ratio_O * 3)

print("--- Analysis of Substance B ---")
print(f"Molar ratios (C:H:O): {ratio_C:.2f} : {ratio_H:.2f} : {ratio_O:.2f}")
print(f"Empirical formula of B: C{emp_C}H{emp_H}O{emp_O}")

# As deduced in the reasoning, structure of B is Pentane-1,3,5-triol
# Its MW is C5H12O3 = 5*12.011 + 12*1.008 + 3*15.999 = 120.12 g/mol
# It has two types of OH groups (primary at C1/C5, secondary at C3), giving two HBr products.
# Product from C3 sub: 3-bromo-1,5-pentanediol -> achiral
# Product from C1/C5 sub: 5-bromo-1,3-pentanediol -> chiral
# This confirms the structure of B.

# Part 2: Analysis of X
# Given combustion data
mass_CO2 = 0.7472  # g
mass_H2O = 0.1834  # g

# Molar masses
M_CO2 = 44.01
M_H2O = 18.015

# Step 1: Calculate moles of C and H in X
moles_C_X = mass_CO2 / M_CO2
moles_H_X = (mass_H2O / M_H2O) * 2

# Step 2: Determine C:H ratio
H_C_ratio = moles_H_X / moles_C_X

# The ratio is ~1.2, which is 6/5. So the formula has a C5H6 unit.
# Possible formulas are (C5H6)nOz

# Step 3: Use molar mass M ~ 150 (range 135-165) to find molecular formula
# n=2: C10H12Oz -> M = 10*12.011 + 12*1.008 + z*15.999 = 132.206 + 16z
z1_mass = 132.206 + 1 * 15.999
z2_mass = 132.206 + 2 * 15.999

print("\n--- Analysis of Substance X ---")
print(f"Moles of C from CO2: {moles_C_X:.4f} mol")
print(f"Moles of H from H2O: {moles_H_X:.4f} mol")
print(f"Ratio of H to C atoms: {H_C_ratio:.2f}")
print("This ratio is approximately 6/5.")
print("With M ~ 150 g/mol, the most plausible molecular formulas are:")
print(f"For C10H12O (z=1), M = {z1_mass:.2f} g/mol")
print(f"For C10H12O2 (z=2), M = {z2_mass:.2f} g/mol")
print("Given the chemical properties (phenol-like and aldehyde-like), C10H12O2 is the more likely candidate.")

final_structure_conclusion = "2-hydroxy-5-propylbenzaldehyde"

print(f"\nFinal Deduced Structure of X is likely an isomer of propyl-hydroxybenzaldehyde, such as {final_structure_conclusion} (formula C10H12O2), as this fits the chemical tests for a phenolic aldehyde and the molecular formula derived from combustion and molar mass data. The ozonolysis information is noted to be inconsistent with simple reaction mechanisms.")
