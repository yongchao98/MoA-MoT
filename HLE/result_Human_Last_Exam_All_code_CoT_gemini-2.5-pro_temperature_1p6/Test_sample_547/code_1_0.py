import math

# --- Given Constants ---
# Molecular Weights (g/mol)
MW_InP = 146
MW_oleate = 281

# Enthalpies (kJ/mol)
H_diss_bulk_InP = -86  # Enthalpy of dissolution for bulk InP
H_protonation_oleate = 7  # Enthalpy of protonation for oleate

# Data for the largest Quantum Dot (QD)
H_diss_large_QD = 70      # Total measured enthalpy of dissolution per mole of InP
mass_frac_oleate_large = 0.46 # 46%

# Data for the smallest Quantum Dot (QD)
H_diss_small_QD = 120     # Total measured enthalpy of dissolution per mole of InP
mass_frac_oleate_small = 0.52 # 52%

# --- Calculations ---
# We want to find the enthalpy contribution from oleate protonation per mole of InP.
# To do this, we first need the mole ratio of oleate to InP for each QD size.

# --- Largest QD Analysis ---
print("--- Analysis for the Largest Quantum Dot ---")
# Assume a 100g sample for calculation
mass_oleate_large = 100 * mass_frac_oleate_large
mass_InP_large = 100 * (1 - mass_frac_oleate_large)

# Calculate moles
moles_oleate_large = mass_oleate_large / MW_oleate
moles_InP_large = mass_InP_large / MW_InP

# Calculate mole ratio (moles of oleate per mole of InP)
mole_ratio_large = moles_oleate_large / moles_InP_large
print(f"Mole ratio (oleate/InP): {mole_ratio_large:.3f}")

# Calculate enthalpy contribution from oleate protonation (per mole of InP)
H_contrib_oleate_large = mole_ratio_large * H_protonation_oleate
print(f"Total Measured Enthalpy: {H_diss_large_QD} kJ/mol of InP")
print(f"Enthalpy from Oleate Protonation = {mole_ratio_large:.3f} mol_oleate/mol_InP * {H_protonation_oleate} kJ/mol_oleate = {H_contrib_oleate_large:.2f} kJ/mol of InP")
print("-" * 40)

# --- Smallest QD Analysis ---
print("--- Analysis for the Smallest Quantum Dot ---")
# Assume a 100g sample
mass_oleate_small = 100 * mass_frac_oleate_small
mass_InP_small = 100 * (1 - mass_frac_oleate_small)

# Calculate moles
moles_oleate_small = mass_oleate_small / MW_oleate
moles_InP_small = mass_InP_small / MW_InP

# Calculate mole ratio (moles of oleate per mole of InP)
mole_ratio_small = moles_oleate_small / moles_InP_small
print(f"Mole ratio (oleate/InP): {mole_ratio_small:.3f}")

# Calculate enthalpy contribution from oleate protonation (per mole of InP)
H_contrib_oleate_small = mole_ratio_small * H_protonation_oleate
print(f"Total Measured Enthalpy: {H_diss_small_QD} kJ/mol of InP")
print(f"Enthalpy from Oleate Protonation = {mole_ratio_small:.3f} mol_oleate/mol_InP * {H_protonation_oleate} kJ/mol_oleate = {H_contrib_oleate_small:.2f} kJ/mol of InP")
print("-" * 40)

# --- Conclusion from Calculation ---
print("--- Comparison and Conclusion ---")
# Calculate the change in total measured enthalpy
delta_H_total = H_diss_small_QD - H_diss_large_QD
print(f"Total observed change in enthalpy: {H_diss_small_QD} - {H_diss_large_QD} = {delta_H_total:.2f} kJ/mol of InP")

# Calculate the change in enthalpy explained by oleate protonation
delta_H_oleate = H_contrib_oleate_small - H_contrib_oleate_large
print(f"Change in enthalpy from oleate protonation: {H_contrib_oleate_small:.2f} - {H_contrib_oleate_large:.2f} = {delta_H_oleate:.2f} kJ/mol of InP")

# Final statement
percent_explained = (delta_H_oleate / delta_H_total) * 100
print(f"The increased oleate protonation only accounts for {delta_H_oleate:.2f} kJ/mol, which is only {percent_explained:.2f}% of the total observed change of {delta_H_total:.2f} kJ/mol.")
print("Therefore, while the effect in option A is real, it is far too small to be the main explanation.")
print("\nEvaluation points to a different, more significant endothermic process that scales with surface area. Option D, the energy needed to disrupt the packed ligand shell, is the most plausible explanation for this large energy difference.")
