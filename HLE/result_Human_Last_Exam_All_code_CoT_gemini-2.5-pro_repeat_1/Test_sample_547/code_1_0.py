# Constants
mw_inp = 146  # g/mol
mw_oleate = 281  # g/mol
h_protonation_oleate = 7  # kJ/mol

# Observed enthalpies
h_dissolution_large_qd = 70  # kJ/mol of InP
h_dissolution_small_qd = 120  # kJ/mol of InP

# Mass fractions for the largest quantum dot
mass_fraction_oleate_large = 0.46
mass_fraction_inp_large = 1.0 - mass_fraction_oleate_large

# Mass fractions for the smallest quantum dot
mass_fraction_oleate_small = 0.52
mass_fraction_inp_small = 1.0 - mass_fraction_oleate_small

# --- Calculations for the largest quantum dot ---

# Assume a 100g sample to find the molar ratio
moles_oleate_large = (100 * mass_fraction_oleate_large) / mw_oleate
moles_inp_large = (100 * mass_fraction_inp_large) / mw_inp
ratio_large = moles_oleate_large / moles_inp_large

# Calculate the enthalpy contribution from oleate protonation (per mole of InP)
h_contrib_oleate_large = ratio_large * h_protonation_oleate

print("--- For the Largest Quantum Dot ---")
print(f"Molar ratio of oleate to InP: {ratio_large:.4f} mol oleate / mol InP")
print(f"Enthalpy contribution from oleate protonation: {ratio_large:.4f} * {h_protonation_oleate} kJ/mol = {h_contrib_oleate_large:.2f} kJ/mol of InP")
print("-" * 35)

# --- Calculations for the smallest quantum dot ---

# Assume a 100g sample to find the molar ratio
moles_oleate_small = (100 * mass_fraction_oleate_small) / mw_oleate
moles_inp_small = (100 * mass_fraction_inp_small) / mw_inp
ratio_small = moles_oleate_small / moles_inp_small

# Calculate the enthalpy contribution from oleate protonation (per mole of InP)
h_contrib_oleate_small = ratio_small * h_protonation_oleate

print("--- For the Smallest Quantum Dot ---")
print(f"Molar ratio of oleate to InP: {ratio_small:.4f} mol oleate / mol InP")
print(f"Enthalpy contribution from oleate protonation: {ratio_small:.4f} * {h_protonation_oleate} kJ/mol = {h_contrib_oleate_small:.2f} kJ/mol of InP")
print("-" * 35)

# --- Comparison and Conclusion ---

# Calculate the change in enthalpy contribution due to oleate protonation
delta_h_oleate = h_contrib_oleate_small - h_contrib_oleate_large

# Calculate the observed change in total dissolution enthalpy
observed_delta_h = h_dissolution_small_qd - h_dissolution_large_qd

print("--- Analysis of the Change ---")
print(f"The observed change in dissolution enthalpy is: {h_dissolution_small_qd} - {h_dissolution_large_qd} = {observed_delta_h} kJ/mol")
print(f"The calculated change due to oleate protonation is: {h_contrib_oleate_small:.2f} - {h_contrib_oleate_large:.2f} = {delta_h_oleate:.2f} kJ/mol")
print("\nConclusion: The enthalpy change from oleate protonation is only a tiny fraction of the observed change. Therefore, explanation A is not sufficient.")
print("The most logical explanation is that another surface-related endothermic process, which scales with the amount of ligand, is dominant. Choice D describes such a process.")
