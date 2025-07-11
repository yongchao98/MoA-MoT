import sys
# Redirect stdout to a variable to prevent printing intermediate steps in the final output
original_stdout = sys.stdout
sys.stdout = None

# --- Step 1: Define given values ---
mw_inp = 146  # g/mol
mw_oleate = 281  # g/mol

dh_diss_inp_bulk = -86  # kJ/mol
dh_prot_oleate = 7  # kJ/mol

# Largest Quantum Dot (QD) data
mass_frac_oleate_large = 0.46  # 46%
dh_obs_large = 70  # kJ/mol

# Smallest Quantum Dot (QD) data
mass_frac_oleate_small = 0.52  # 52%
dh_obs_small = 120  # kJ/mol

# --- Step 2: Calculate molar ratios for a 100g sample ---

# For the largest QD
mass_oleate_large = 100 * mass_frac_oleate_large
mass_inp_large = 100 - mass_oleate_large
moles_oleate_large = mass_oleate_large / mw_oleate
moles_inp_large = mass_inp_large / mw_inp
ratio_large = moles_oleate_large / moles_inp_large

# For the smallest QD
mass_oleate_small = 100 * mass_frac_oleate_small
mass_inp_small = 100 - mass_oleate_small
moles_oleate_small = mass_oleate_small / mw_oleate
moles_inp_small = mass_inp_small / mw_inp
ratio_small = moles_oleate_small / moles_inp_small

# --- Step 3: Calculate the expected enthalpy from known reactions ---

# Expected enthalpy if only InP dissolution and oleate protonation occur
# Formula: ΔH_expected = ΔH_InP + (ratio_oleate/InP) * ΔH_oleate
dh_expected_large = dh_diss_inp_bulk + (ratio_large * dh_prot_oleate)
dh_expected_small = dh_diss_inp_bulk + (ratio_small * dh_prot_oleate)

# --- Step 4: Calculate the unaccounted-for enthalpy ---

# Unaccounted enthalpy = Observed enthalpy - Expected enthalpy
dh_unaccounted_large = dh_obs_large - dh_expected_large
dh_unaccounted_small = dh_obs_small - dh_expected_small

# --- Step 5: Print the results and analysis ---
sys.stdout = original_stdout # Restore stdout

print("Analysis of Enthalpy Contributions (all values in kJ/mol of InP)\n")

print("--- For the Largest Quantum Dot ---")
print(f"The molar ratio of oleate to InP is: {ratio_large:.2f}")
print(f"The enthalpy contribution from oleate protonation is: {ratio_large:.2f} * {dh_prot_oleate} kJ/mol = {ratio_large * dh_prot_oleate:.2f} kJ/mol")
print(f"The total expected enthalpy from known reactions (InP dissolution + oleate protonation) is: {dh_diss_inp_bulk} + {ratio_large * dh_prot_oleate:.2f} = {dh_expected_large:.2f} kJ/mol")
print(f"The observed enthalpy is: {dh_obs_large} kJ/mol")
print(f"The unaccounted-for endothermic enthalpy is: {dh_obs_large} - ({dh_expected_large:.2f}) = {dh_unaccounted_large:.2f} kJ/mol\n")


print("--- For the Smallest Quantum Dot ---")
print(f"The molar ratio of oleate to InP is: {ratio_small:.2f}")
print(f"The enthalpy contribution from oleate protonation is: {ratio_small:.2f} * {dh_prot_oleate} kJ/mol = {ratio_small * dh_prot_oleate:.2f} kJ/mol")
print(f"The total expected enthalpy from known reactions (InP dissolution + oleate protonation) is: {dh_diss_inp_bulk} + {ratio_small * dh_prot_oleate:.2f} = {dh_expected_small:.2f} kJ/mol")
print(f"The observed enthalpy is: {dh_obs_small} kJ/mol")
print(f"The unaccounted-for endothermic enthalpy is: {dh_obs_small} - ({dh_expected_small:.2f}) = {dh_unaccounted_small:.2f} kJ/mol\n")

print("--- Conclusion ---")
print("The analysis shows that the protonation of oleate (Choice A) only contributes a small amount to the total enthalpy ({:.2f} to {:.2f} kJ/mol).".format(ratio_large * dh_prot_oleate, ratio_small * dh_prot_oleate))
print("There is a large, unexplained endothermic contribution (+{:.2f} kJ/mol for large dots, +{:.2f} kJ/mol for small dots).".format(dh_unaccounted_large, dh_unaccounted_small))
print("This unexplained energy is larger for the smaller quantum dots, which have a higher surface area and a greater proportion of ligands.")
print("This strongly suggests that a significant amount of energy is required to break the interactions within the packed organic ligand shell on the surface, which is the explanation provided in Choice D.")