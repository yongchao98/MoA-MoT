# Molar masses
M_C = 12.011
M_H = 1.008
M_CO2 = 44.01
M_H2O = 18.015

# Given masses of products
mass_CO2 = 0.7472  # g
mass_H2O = 0.1834  # g

# --- Step 1: Calculate moles of C and H ---
# Moles of Carbon from CO2
moles_C = mass_CO2 / M_CO2
# Moles of Hydrogen from H2O (note the factor of 2)
moles_H = (mass_H2O / M_H2O) * 2

# --- Step 2: Determine the simplest ratio ---
# To find the simplest integer ratio, divide both by the smaller number
ratio_H_to_C = moles_H / moles_C

# To express this as a simple fraction (like 6/5)
# We can test common fractions. 6/5 = 1.2
# The calculated ratio is very close to 1.2, so the empirical ratio is C:H = 5:6

print(f"Calculation of Empirical Ratio C:H")
print(f"Mass of CO2 formed = {mass_CO2} g")
print(f"Mass of H2O formed = {mass_H2O} g")
print(f"Moles of Carbon (C) = {moles_C:.4f} mol")
print(f"Moles of Hydrogen (H) = {moles_H:.4f} mol")
print(f"Ratio of H moles to C moles = {ratio_H_to_C:.4f}")
print(f"This corresponds to an empirical formula ratio for C:H of 5:6.")
print("\nFinal Proposed Structure for X (based on analysis of all data):")
print("Structure Name: 3,4-diacetylhex-3-ene-2,5-dione")
print("Molecular Formula: C10H12O4")
# The final equation is the structure name and formula
print("Equation: X = 3,4-diacetylhex-3-ene-2,5-dione")
print("C = 10")
print("H = 12")
print("O = 4")