import sys
# Redirect stdout to a string buffer to capture all prints
from io import StringIO
old_stdout = sys.stdout
sys.stdout = mystdout = StringIO()

# --- GIVEN CONSTANTS ---
MW_InP = 146  # g/mol
MW_oleate = 281  # g/mol
H_prot_oleate = 7  # kJ/mol of oleate
H_diss_large_QD_total = 70  # kJ/mol of InP
H_diss_small_QD_total = 120  # kJ/mol of InP
mass_frac_oleate_large = 0.46
mass_frac_oleate_small = 0.52

# --- CALCULATIONS FOR THE LARGEST QUANTUM DOT ---

# Based on a 100g sample, calculate the moles of each component
mass_InP_large = 100 * (1 - mass_frac_oleate_large)
mass_oleate_large = 100 * mass_frac_oleate_large
moles_InP_large = mass_InP_large / MW_InP
moles_oleate_large = mass_oleate_large / MW_oleate

# Determine the molar ratio of oleate to InP
ratio_large = moles_oleate_large / moles_InP_large

# Calculate the enthalpy contribution from oleate protonation (per mole of InP)
H_contrib_oleate_large = ratio_large * H_prot_oleate


# --- CALCULATIONS FOR THE SMALLEST QUANTUM DOT ---

# Based on a 100g sample, calculate the moles of each component
mass_InP_small = 100 * (1 - mass_frac_oleate_small)
mass_oleate_small = 100 * mass_frac_oleate_small
moles_InP_small = mass_InP_small / MW_InP
moles_oleate_small = mass_oleate_small / MW_oleate

# Determine the molar ratio of oleate to InP
ratio_small = moles_oleate_small / moles_InP_small

# Calculate the enthalpy contribution from oleate protonation (per mole of InP)
H_contrib_oleate_small = ratio_small * H_prot_oleate


# --- ANALYSIS OF CHANGES ---

# Calculate the total observed change in dissolution enthalpy
total_enthalpy_change = H_diss_small_QD_total - H_diss_large_QD_total

# Calculate the change in enthalpy attributable solely to oleate protonation
oleate_enthalpy_change = H_contrib_oleate_small - H_contrib_oleate_large

# --- PRINTING THE RESULTS ---

print("Quantitative Analysis of Enthalpy Contributions")
print("=================================================")

print("\n--- Largest Quantum Dot ---")
print(f"The molar ratio of oleate to InP is {moles_oleate_large:.3f} / {moles_InP_large:.3f} = {ratio_large:.3f} mol/mol.")
print(f"The endothermic contribution from oleate protonation is {ratio_large:.3f} * {H_prot_oleate} kJ/mol = {H_contrib_oleate_large:.2f} kJ per mole of InP.")

print("\n--- Smallest Quantum Dot ---")
print(f"The molar ratio of oleate to InP is {moles_oleate_small:.3f} / {moles_InP_small:.3f} = {ratio_small:.3f} mol/mol.")
print(f"The endothermic contribution from oleate protonation is {ratio_small:.3f} * {H_prot_oleate} kJ/mol = {H_contrib_oleate_small:.2f} kJ per mole of InP.")

print("\n--- Comparison of Changes ---")
print(f"The total measured enthalpy of dissolution increases from {H_diss_large_QD_total} kJ/mol to {H_diss_small_QD_total} kJ/mol.")
print(f"Total observed increase in endothermicity = {H_diss_small_QD_total} - {H_diss_large_QD_total} = {total_enthalpy_change:.2f} kJ/mol.")
print(f"Increase in endothermicity from oleate protonation = {H_contrib_oleate_small:.2f} - {H_contrib_oleate_large:.2f} = {oleate_enthalpy_change:.2f} kJ/mol.")

print("\n--- Conclusion ---")
print(f"The analysis shows that the increased amount of oleate in smaller dots makes the dissolution more endothermic by only {oleate_enthalpy_change:.2f} kJ/mol.")
print(f"This accounts for only a tiny fraction ({oleate_enthalpy_change/total_enthalpy_change:.1%}) of the total observed increase of {total_enthalpy_change:.2f} kJ/mol.")
print("Therefore, the protonation of oleate (Option A) is not the most logical explanation for the large observed change.")
print("A much larger endothermic effect that scales with the amount of ligand is required, such as disrupting the interactions within the tightly packed ligand shell (Option D).")

# Restore original stdout and print the captured output
sys.stdout = old_stdout
print(mystdout.getvalue())