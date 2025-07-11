import math

# --- Given Constants ---
MW_InP = 146  # g/mol
MW_oleate = 281  # g/mol
dH_protonation_oleate = 7  # kJ/mol of oleate
dH_diss_bulk_InP = -86 # kJ/mol of InP

# --- Data for the Largest Quantum Dot ---
mass_frac_oleate_large = 0.46  # 46%
dH_diss_large_QD = 70  # kJ/mol of InP

# --- Data for the Smallest Quantum Dot ---
mass_frac_oleate_small = 0.52  # 52%
dH_diss_small_QD = 120  # kJ/mol of InP

# --- Calculations ---

# Assume a 100g sample for calculations
sample_mass = 100  # g

# For the Largest QD
mass_oleate_large = sample_mass * mass_frac_oleate_large
mass_InP_large = sample_mass * (1 - mass_frac_oleate_large)
moles_oleate_large = mass_oleate_large / MW_oleate
moles_InP_large = mass_InP_large / MW_InP
mole_ratio_large = moles_oleate_large / moles_InP_large
dH_contribution_large = mole_ratio_large * dH_protonation_oleate

# For the Smallest QD
mass_oleate_small = sample_mass * mass_frac_oleate_small
mass_InP_small = sample_mass * (1 - mass_frac_oleate_small)
moles_oleate_small = mass_oleate_small / MW_oleate
moles_InP_small = mass_InP_small / MW_InP
mole_ratio_small = moles_oleate_small / moles_InP_small
dH_contribution_small = mole_ratio_small * dH_protonation_oleate

# --- Comparison ---
observed_dH_difference = dH_diss_small_QD - dH_diss_large_QD
calculated_dH_difference_from_oleate = dH_contribution_small - dH_contribution_large

# --- Output the analysis ---
print("--- Analysis of Enthalpy Contributions ---")
print("\n--- For the Largest Quantum Dot ---")
print(f"For every 1 mole of InP, there are {mole_ratio_large:.2f} moles of oleate.")
print(f"The enthalpy contribution from oleate protonation is: {mole_ratio_large:.2f} * {dH_protonation_oleate} kJ/mol = {dH_contribution_large:.2f} kJ per mole of InP.")
print(f"Equation: {dH_diss_large_QD} kJ/mol (total) = [InP dissolution] + {dH_contribution_large:.2f} kJ/mol (oleate protonation)")


print("\n--- For the Smallest Quantum Dot ---")
print(f"For every 1 mole of InP, there are {mole_ratio_small:.2f} moles of oleate.")
print(f"The enthalpy contribution from oleate protonation is: {mole_ratio_small:.2f} * {dH_protonation_oleate} kJ/mol = {dH_contribution_small:.2f} kJ per mole of InP.")
print(f"Equation: {dH_diss_small_QD} kJ/mol (total) = [InP dissolution] + {dH_contribution_small:.2f} kJ/mol (oleate protonation)")


print("\n--- Comparison and Conclusion ---")
print(f"The observed increase in dissolution enthalpy is: {dH_diss_small_QD} - {dH_diss_large_QD} = {observed_dH_difference:.2f} kJ/mol of InP.")
print(f"The calculated increase from oleate protonation alone is: {dH_contribution_small:.2f} - {dH_contribution_large:.2f} = {calculated_dH_difference_from_oleate:.2f} kJ/mol of InP.")
print("\nConclusion: The enthalpy change due to the increased amount of oleate ({:.2f} kJ/mol) is far too small to explain the total observed enthalpy change ({:.2f} kJ/mol).".format(calculated_dH_difference_from_oleate, observed_dH_difference))
print("This suggests another, more significant endothermic process that scales with surface area must be responsible. The disruption of the packed ligand shell (Choice D) is the most plausible explanation for this large energy difference.")
