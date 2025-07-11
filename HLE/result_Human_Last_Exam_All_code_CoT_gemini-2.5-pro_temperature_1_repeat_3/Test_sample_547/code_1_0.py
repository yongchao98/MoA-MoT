# Given values
mw_inp = 146  # g/mol
mw_oleate = 281  # g/mol
h_protonation_oleate = 7  # kJ/mol

# Largest quantum dot data
mass_fraction_oleate_large = 0.46
h_diss_large = 70  # kJ/mol of InP

# Smallest quantum dot data
mass_fraction_oleate_small = 0.52
h_diss_small = 120  # kJ/mol of InP

# --- Calculations for the largest quantum dot ---
# Assume a 100g sample for calculation
mass_oleate_large = 100 * mass_fraction_oleate_large
mass_inp_large = 100 * (1 - mass_fraction_oleate_large)

# Moles of each component in the 100g sample
moles_oleate_large = mass_oleate_large / mw_oleate
moles_inp_large = mass_inp_large / mw_inp

# Molar ratio of oleate to InP
ratio_large = moles_oleate_large / moles_inp_large

# Enthalpy contribution from oleate protonation per mole of InP
h_from_oleate_large = ratio_large * h_protonation_oleate

# --- Calculations for the smallest quantum dot ---
# Assume a 100g sample for calculation
mass_oleate_small = 100 * mass_fraction_oleate_small
mass_inp_small = 100 * (1 - mass_fraction_oleate_small)

# Moles of each component in the 100g sample
moles_oleate_small = mass_oleate_small / mw_oleate
moles_inp_small = mass_inp_small / mw_inp

# Molar ratio of oleate to InP
ratio_small = moles_oleate_small / moles_inp_small

# Enthalpy contribution from oleate protonation per mole of InP
h_from_oleate_small = ratio_small * h_protonation_oleate

# --- Analysis of the change in enthalpy ---
# Observed change in dissolution enthalpy
observed_h_change = h_diss_small - h_diss_large

# Calculated change in enthalpy due to oleate protonation
calculated_h_change_from_oleate = h_from_oleate_small - h_from_oleate_large

# --- Output the results ---
print("Analysis of Enthalpy Contributions:\n")

print(f"For the largest quantum dot:")
print(f"  - The molar ratio of oleate to InP is: {ratio_large:.3f}")
print(f"  - The enthalpy from oleate protonation is {ratio_large:.3f} * {h_protonation_oleate} kJ/mol = {h_from_oleate_large:.2f} kJ/mol of InP")
print("-" * 30)

print(f"For the smallest quantum dot:")
print(f"  - The molar ratio of oleate to InP is: {ratio_small:.3f}")
print(f"  - The enthalpy from oleate protonation is {ratio_small:.3f} * {h_protonation_oleate} kJ/mol = {h_from_oleate_small:.2f} kJ/mol of InP")
print("-" * 30)

print("Comparison:")
print(f"The total observed change in dissolution enthalpy is: {h_diss_small} - {h_diss_large} = {observed_h_change} kJ/mol")
print(f"The change in enthalpy from oleate protonation is: {h_from_oleate_small:.2f} - {h_from_oleate_large:.2f} = {calculated_h_change_from_oleate:.2f} kJ/mol")
print("\nConclusion:")
print(f"The contribution from oleate protonation ({calculated_h_change_from_oleate:.2f} kJ/mol) is only a small fraction of the total observed change ({observed_h_change} kJ/mol).")
print("Therefore, the increased amount of oleate and its protonation (Option A) cannot be the main explanation for the large increase in endothermicity.")
