# Constants from the problem description
mw_inp = 146  # g/mol
mw_oleate = 281  # g/mol
h_protonation_oleate = 7  # kJ/mol

# Data for the largest quantum dot (QD)
h_diss_large = 70  # kJ/mol of InP
mass_frac_oleate_large = 0.46
mass_frac_inp_large = 1 - mass_frac_oleate_large

# Data for the smallest quantum dot (QD)
h_diss_small = 120  # kJ/mol of InP
mass_frac_oleate_small = 0.52
mass_frac_inp_small = 1 - mass_frac_oleate_small

# --- Calculations for the largest QD ---
# Assume a 100g sample to find the molar ratio
moles_oleate_large = (100 * mass_frac_oleate_large) / mw_oleate
moles_inp_large = (100 * mass_frac_inp_large) / mw_inp
molar_ratio_large = moles_oleate_large / moles_inp_large
h_contrib_protonation_large = molar_ratio_large * h_protonation_oleate

# --- Calculations for the smallest QD ---
# Assume a 100g sample to find the molar ratio
moles_oleate_small = (100 * mass_frac_oleate_small) / mw_oleate
moles_inp_small = (100 * mass_frac_inp_small) / mw_inp
molar_ratio_small = moles_oleate_small / moles_inp_small
h_contrib_protonation_small = molar_ratio_small * h_protonation_oleate

# --- Comparison ---
# Total observed change in enthalpy
total_enthalpy_change = h_diss_small - h_diss_large

# Change in enthalpy due to oleate protonation
protonation_enthalpy_change = h_contrib_protonation_small - h_contrib_protonation_large

# --- Output the results ---
print("Analysis of Enthalpy Contributions:")
print(f"Enthalpy contribution from oleate protonation for largest QD: {h_contrib_protonation_large:.2f} kJ/mol of InP")
print(f"Enthalpy contribution from oleate protonation for smallest QD: {h_contrib_protonation_small:.2f} kJ/mol of InP")
print("-" * 30)
print(f"Total observed change in dissolution enthalpy: {h_diss_small} kJ/mol - {h_diss_large} kJ/mol = {total_enthalpy_change:.2f} kJ/mol")
print(f"Change in enthalpy explained by protonation: {h_contrib_protonation_small:.2f} kJ/mol - {h_contrib_protonation_large:.2f} kJ/mol = {protonation_enthalpy_change:.2f} kJ/mol")
print("-" * 30)
print("Conclusion:")
print(f"The change due to protonation ({protonation_enthalpy_change:.2f} kJ/mol) is significantly smaller than the total observed change ({total_enthalpy_change:.2f} kJ/mol).")
print("Therefore, the increased amount of oleate being protonated (Answer A) is not the main cause of the large endothermic shift.")
print("This suggests another, more significant endothermic process is at play, which scales with surface area. Disrupting the tightly packed ligand shell (Answer D) is the most plausible explanation.")
