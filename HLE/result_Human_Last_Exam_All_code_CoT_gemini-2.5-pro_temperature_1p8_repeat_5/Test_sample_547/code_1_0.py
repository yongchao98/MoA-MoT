import sys
import io

# Backup the original stdout
original_stdout = sys.stdout
# Create an in-memory stream for capturing the output
captured_output = io.StringIO()
# Redirect stdout to the in-memory stream
sys.stdout = captured_output

# --- Given Constants ---
# Molecular weights (g/mol)
mw_inp = 146
mw_oleate = 281

# Enthalpies (kJ/mol)
dh_protonation_oleate = 7  # Per mole of oleate
dh_diss_large_qd = 70      # Per mole of InP
dh_diss_small_qd = 120     # Per mole of InP

# Mass fractions
mass_fraction_oleate_large = 0.46
mass_fraction_oleate_small = 0.52

# --- Calculations for the Largest Quantum Dot ---
# For a 100g sample
mass_oleate_large = 100 * mass_fraction_oleate_large
mass_inp_large = 100 - mass_oleate_large

# Moles
moles_oleate_large = mass_oleate_large / mw_oleate
moles_inp_large = mass_inp_large / mw_inp

# Molar ratio of oleate to InP
molar_ratio_large = moles_oleate_large / moles_inp_large

# Enthalpy contribution from oleate protonation (per mole of InP)
dh_protonation_contribution_large = molar_ratio_large * dh_protonation_oleate


# --- Calculations for the Smallest Quantum Dot ---
# For a 100g sample
mass_oleate_small = 100 * mass_fraction_oleate_small
mass_inp_small = 100 - mass_oleate_small

# Moles
moles_oleate_small = mass_oleate_small / mw_oleate
moles_inp_small = mass_inp_small / mw_inp

# Molar ratio of oleate to InP
molar_ratio_small = moles_oleate_small / moles_inp_small

# Enthalpy contribution from oleate protonation (per mole of InP)
dh_protonation_contribution_small = molar_ratio_small * dh_protonation_oleate

# --- Comparison and Conclusion ---
# Observed total change in dissolution enthalpy
total_dh_change_observed = dh_diss_small_qd - dh_diss_large_qd

# Calculated change in enthalpy due to oleate protonation
dh_change_from_protonation = dh_protonation_contribution_small - dh_protonation_contribution_large

# --- Output the results ---
print("Analysis of Enthalpy Contributions:")
print("-" * 40)
print(f"Contribution from oleate protonation (Largest QD): {dh_protonation_contribution_large:.2f} kJ/mol of InP")
print(f"Contribution from oleate protonation (Smallest QD): {dh_protonation_contribution_small:.2f} kJ/mol of InP")
print("-" * 40)

print("Observed total change in dissolution enthalpy:")
print(f"ΔH_small - ΔH_large = {dh_diss_small_qd} kJ/mol - {dh_diss_large_qd} kJ/mol = {total_dh_change_observed:.2f} kJ/mol")
print("\nCalculated enthalpy change due to oleate protonation:")
print(f"Contribution_small - Contribution_large = {dh_protonation_contribution_small:.2f} kJ/mol - {dh_protonation_contribution_large:.2f} kJ/mol = {dh_change_from_protonation:.2f} kJ/mol")

print("\nConclusion:")
print(f"The change in oleate protonation enthalpy ({dh_change_from_protonation:.2f} kJ/mol) accounts for only a small fraction of the total observed change ({total_dh_change_observed:.2f} kJ/mol).")
print("Therefore, oleate protonation is not the primary explanation for the trend.")

# --- Finalizing the output ---
# Restore the original stdout
sys.stdout = original_stdout
# Get the captured output as a string
output_string = captured_output.getvalue()
# Print the captured output
print(output_string)

# Analysis based on the calculation:
# The calculation shows that the change in enthalpy from protonating the different amounts of oleate is only ~0.83 kJ/mol, whereas the observed change is 50 kJ/mol. This disproves answer A.
#
# Let's briefly analyze other options:
# B. The enthalpies are normalized per mole of InP, so the amount of InP is irrelevant. Incorrect.
# C. A positive surface energy of the InP core would make the dissolution *more* exothermic (less endothermic), the opposite of the observation. Incorrect.
# E. This is speculative and less direct than considering the main structural difference: the ligand shell.
# D. The dissolution requires breaking the bonds between the InP surface and the oleate ligands, as well as the van der Waals interactions between the packed oleate tails. This is an endothermic process. Smaller QDs have a higher surface-to-volume ratio, meaning a greater proportion of the system's energy is stored in this surface/ligand interface. Disrupting this interface becomes a larger energetic cost per mole of InP for smaller dots. This aligns perfectly with the observation of a large and size-dependent endothermic enthalpy. This is the most logical explanation.