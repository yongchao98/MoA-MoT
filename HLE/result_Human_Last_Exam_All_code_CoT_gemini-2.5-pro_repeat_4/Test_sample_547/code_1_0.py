import sys

# --- Given Constants ---
# Molecular Weights (g/mol)
mw_inp = 146
mw_oleate = 281

# Enthalpies (kJ/mol)
# The enthalpy of dissolution for bulk InP is -86 kJ/mol.
# The measured enthalpies for the QDs are given per mole of InP.
h_diss_large_qd = 70  # Enthalpy of dissolution for the largest QD
h_diss_small_qd = 120 # Enthalpy of dissolution for the smallest QD
h_prot_oleate = 7     # Enthalpy of protonation for oleate

# Mass Fractions
mass_frac_oleate_large = 0.46 # 46% for the largest QD
mass_frac_oleate_small = 0.52 # 52% for the smallest QD

# --- Plan ---
# 1. Calculate the mole ratio of oleate to InP for each QD size.
#    We assume a 100g sample for ease of calculation.

# For the largest QD
mass_inp_large = 100 * (1 - mass_frac_oleate_large)
mass_oleate_large = 100 * mass_frac_oleate_large
moles_inp_large = mass_inp_large / mw_inp
moles_oleate_large = mass_oleate_large / mw_oleate
mole_ratio_large = moles_oleate_large / moles_inp_large

# For the smallest QD
mass_inp_small = 100 * (1 - mass_frac_oleate_small)
mass_oleate_small = 100 * mass_frac_oleate_small
moles_inp_small = mass_inp_small / mw_inp
moles_oleate_small = mass_oleate_small / mw_oleate
mole_ratio_small = moles_oleate_small / moles_inp_small

# 2. Calculate the enthalpy contribution from oleate protonation (per mole of InP).
h_contrib_prot_large = mole_ratio_large * h_prot_oleate
h_contrib_prot_small = mole_ratio_small * h_prot_oleate

# 3. Compare the change in protonation enthalpy to the total observed change.
total_h_change_observed = h_diss_small_qd - h_diss_large_qd
h_change_from_protonation = h_contrib_prot_small - h_contrib_prot_large

# --- Output the results and reasoning ---
print("Analysis of Enthalpy Contributions:\n")

print(f"For the largest quantum dot:")
print(f"  - Moles of Oleate per mole of InP = ({mass_frac_oleate_large} / {mw_oleate}) / ((1 - {mass_frac_oleate_large}) / {mw_inp}) = {mole_ratio_large:.4f}")
print(f"  - Enthalpy from Oleate Protonation = {mole_ratio_large:.4f} mol * {h_prot_oleate} kJ/mol = {h_contrib_prot_large:.2f} kJ per mole of InP\n")

print(f"For the smallest quantum dot:")
print(f"  - Moles of Oleate per mole of InP = ({mass_frac_oleate_small} / {mw_oleate}) / ((1 - {mass_frac_oleate_small}) / {mw_inp}) = {mole_ratio_small:.4f}")
print(f"  - Enthalpy from Oleate Protonation = {mole_ratio_small:.4f} mol * {h_prot_oleate} kJ/mol = {h_contrib_prot_small:.2f} kJ per mole of InP\n")

print("Comparison of Enthalpy Changes:")
print("The final equation compares the observed change to the change explained by oleate protonation.\n")
print(f"Total Observed Change in Enthalpy = {h_diss_small_qd} kJ/mol - {h_diss_large_qd} kJ/mol = {total_h_change_observed} kJ/mol")
print(f"Change Explained by Oleate Protonation = {h_contrib_prot_small:.2f} kJ/mol - {h_contrib_prot_large:.2f} kJ/mol = {h_change_from_protonation:.2f} kJ/mol\n")

print("Conclusion:")
print(f"The increasing amount of oleate in smaller quantum dots makes the dissolution more endothermic, but only by {h_change_from_protonation:.2f} kJ/mol.")
print(f"This is a very small fraction of the total observed endothermic shift of {total_h_change_observed} kJ/mol.")
print("Therefore, while oleate protonation contributes, it cannot be the primary explanation (eliminating Answer A).")
print("The most logical explanation is D, which posits that a large amount of energy is required to disrupt the tightly packed shell of organic ligands on the QD surface, an effect that becomes larger as the surface-to-volume ratio increases for smaller dots.")
