# Constants and given data
MW_InP = 146  # g/mol
MW_oleate = 281  # g/mol
dH_protonation_oleate = 7  # kJ/mol of oleate

# Data for the largest quantum dot (LQD)
mass_frac_oleate_LQD = 0.46
dH_dissolution_LQD = 70  # kJ/mol of InP

# Data for the smallest quantum dot (SQD)
mass_frac_oleate_SQD = 0.52
dH_dissolution_SQD = 120  # kJ/mol of InP

# --- Step 1: Analyze the Largest Quantum Dot ---
# In a hypothetical 100g sample of LQD material:
mass_oleate_LQD = mass_frac_oleate_LQD * 100
mass_InP_LQD = (1 - mass_frac_oleate_LQD) * 100
moles_oleate_LQD = mass_oleate_LQD / MW_oleate
moles_InP_LQD = mass_InP_LQD / MW_InP
# Calculate the molar ratio of oleate to InP
ratio_LQD = moles_oleate_LQD / moles_InP_LQD
# Calculate the enthalpy contribution from oleate protonation per mole of InP
dH_from_oleate_LQD = ratio_LQD * dH_protonation_oleate

# --- Step 2: Analyze the Smallest Quantum Dot ---
# In a hypothetical 100g sample of SQD material:
mass_oleate_SQD = mass_frac_oleate_SQD * 100
mass_InP_SQD = (1 - mass_frac_oleate_SQD) * 100
moles_oleate_SQD = mass_oleate_SQD / MW_oleate
moles_InP_SQD = mass_InP_SQD / MW_InP
# Calculate the molar ratio of oleate to InP
ratio_SQD = moles_oleate_SQD / moles_InP_SQD
# Calculate the enthalpy contribution from oleate protonation per mole of InP
dH_from_oleate_SQD = ratio_SQD * dH_protonation_oleate

# --- Step 3: Compare the calculated effect with the observed effect ---
# The total observed change in dissolution enthalpy
observed_dH_change = dH_dissolution_SQD - dH_dissolution_LQD

# The change in enthalpy that can be explained by oleate protonation alone
calculated_dH_change_from_oleate = dH_from_oleate_SQD - dH_from_oleate_LQD

# --- Step 4: Print the results and analysis ---
print("This analysis tests the hypothesis that the increase in dissolution enthalpy is caused by the protonation of oleate ligands (Choice A).")
print("-" * 80)

print(f"For the largest quantum dot:")
print(f"  The molar ratio of oleate to InP is {ratio_LQD:.3f} mol/mol.")
print(f"  The endothermic contribution from oleate protonation is {ratio_LQD:.3f} * {dH_protonation_oleate} kJ/mol = {dH_from_oleate_LQD:.2f} kJ/mol of InP.")
print()

print(f"For the smallest quantum dot:")
print(f"  The molar ratio of oleate to InP is {ratio_SQD:.3f} mol/mol.")
print(f"  The endothermic contribution from oleate protonation is {ratio_SQD:.3f} * {dH_protonation_oleate} kJ/mol = {dH_from_oleate_SQD:.2f} kJ/mol of InP.")
print("-" * 80)

print("Final Comparison:")
print(f"The observed total increase in enthalpy of dissolution is {dH_dissolution_SQD} - {dH_dissolution_LQD} = {observed_dH_change:.2f} kJ/mol.")
print(f"The calculated increase in enthalpy from oleate protonation alone is {dH_from_oleate_SQD:.2f} - {dH_from_oleate_LQD:.2f} = {calculated_dH_change_from_oleate:.2f} kJ/mol.")
print()
print(f"Conclusion: The additional oleate protonation only accounts for {calculated_dH_change_from_oleate:.2f} kJ/mol of the observed {observed_dH_change:.2f} kJ/mol change. This is only {calculated_dH_change_from_oleate/observed_dH_change:.1%} of the total effect.")
print("Therefore, another, much larger endothermic process that scales with the amount of ligand must be responsible. Choice D provides this explanation.")

<<<D>>>