import sys
import io

# Backup the original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

# --- Given Constants ---
# Molecular Weights (g/mol)
MW_InP = 146
MW_oleate = 281

# Enthalpies (kJ/mol)
H_diss_bulk_InP = -86
H_protonation_oleate = 7
H_diss_largest_QD = 70
H_diss_smallest_QD = 120

# Mass Fractions (%)
mass_frac_oleate_largest = 0.46
mass_frac_oleate_smallest = 0.52

# --- Calculations for the Largest Quantum Dot ---
print("--- Analysis for the Largest Quantum Dot ---")
# In a 100g sample:
mass_oleate_largest = 100 * mass_frac_oleate_largest
mass_InP_largest = 100 * (1 - mass_frac_oleate_largest)

# Moles in the 100g sample
moles_oleate_largest = mass_oleate_largest / MW_oleate
moles_InP_largest = mass_InP_largest / MW_InP

# Molar ratio of oleate to InP
ratio_largest = moles_oleate_largest / moles_InP_largest
print(f"Molar ratio (oleate/InP) for largest QD: {ratio_largest:.3f}")

# Enthalpy contribution from oleate protonation (per mole of InP)
H_oleate_contrib_largest = ratio_largest * H_protonation_oleate
print(f"Enthalpy from oleate protonation (per mole InP): {ratio_largest:.3f} mol_oleate/mol_InP * {H_protonation_oleate} kJ/mol_oleate = {H_oleate_contrib_largest:.2f} kJ/mol_InP")
print("")

# --- Calculations for the Smallest Quantum Dot ---
print("--- Analysis for the Smallest Quantum Dot ---")
# In a 100g sample:
mass_oleate_smallest = 100 * mass_frac_oleate_smallest
mass_InP_smallest = 100 * (1 - mass_frac_oleate_smallest)

# Moles in the 100g sample
moles_oleate_smallest = mass_oleate_smallest / MW_oleate
moles_InP_smallest = mass_InP_smallest / MW_InP

# Molar ratio of oleate to InP
ratio_smallest = moles_oleate_smallest / moles_InP_smallest
print(f"Molar ratio (oleate/InP) for smallest QD: {ratio_smallest:.3f}")

# Enthalpy contribution from oleate protonation (per mole of InP)
H_oleate_contrib_smallest = ratio_smallest * H_protonation_oleate
print(f"Enthalpy from oleate protonation (per mole InP): {ratio_smallest:.3f} mol_oleate/mol_InP * {H_protonation_oleate} kJ/mol_oleate = {H_oleate_contrib_smallest:.2f} kJ/mol_InP")
print("")

# --- Comparison and Conclusion ---
print("--- Comparison of Changes ---")
# Total change in measured dissolution enthalpy
total_H_change = H_diss_smallest_QD - H_diss_largest_QD
print(f"Total observed change in enthalpy: {H_diss_smallest_QD} kJ/mol - {H_diss_largest_QD} kJ/mol = {total_H_change} kJ/mol")

# Change in enthalpy due to oleate protonation
H_change_from_oleate = H_oleate_contrib_smallest - H_oleate_contrib_largest
print(f"Change in enthalpy from oleate protonation: {H_oleate_contrib_smallest:.2f} kJ/mol - {H_oleate_contrib_largest:.2f} kJ/mol = {H_change_from_oleate:.2f} kJ/mol")
print("")
print("Conclusion:")
print(f"The protonation of the additional oleate on the smaller dots only accounts for {H_change_from_oleate:.2f} kJ/mol of the total {total_H_change} kJ/mol increase in endothermicity.")
print("This is a very small fraction of the total effect. Therefore, another, more significant endothermic process that scales with the surface area must be responsible.")
print("This points to the energy required to disrupt the tightly packed shell of organic ligands on the quantum dot surface, as described in choice D.")

# --- Final Answer ---
# Restore the original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()
# Print the captured output to the actual console
print(output)