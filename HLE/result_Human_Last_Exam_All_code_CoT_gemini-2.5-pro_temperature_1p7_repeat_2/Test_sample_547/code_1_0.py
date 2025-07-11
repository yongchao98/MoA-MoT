import sys

# Define constants based on the problem description
MW_InP = 146.0  # g/mol
MW_oleate = 281.0  # g/mol
dH_prot_oleate = 7.0  # kJ/mol

# --- Data for the largest quantum dot ---
mass_frac_oleate_large = 0.46
mass_frac_InP_large = 1.0 - mass_frac_oleate_large
dH_diss_large_measured = 70.0  # kJ/mol of InP

# --- Data for the smallest quantum dot ---
mass_frac_oleate_small = 0.52
mass_frac_InP_small = 1.0 - mass_frac_oleate_small
dH_diss_small_measured = 120.0  # kJ/mol of InP

# --- Calculation for the largest QD ---
# Molar ratio of oleate to InP in the largest QD
# Consider a 100g sample
moles_oleate_large = (mass_frac_oleate_large * 100) / MW_oleate
moles_InP_large = (mass_frac_InP_large * 100) / MW_InP
molar_ratio_large = moles_oleate_large / moles_InP_large

# Enthalpy contribution from oleate protonation (per mole of InP)
dH_from_oleate_large = molar_ratio_large * dH_prot_oleate

# --- Calculation for the smallest QD ---
# Molar ratio of oleate to InP in the smallest QD
# Consider a 100g sample
moles_oleate_small = (mass_frac_oleate_small * 100) / MW_oleate
moles_InP_small = (mass_frac_InP_small * 100) / MW_InP
molar_ratio_small = moles_oleate_small / moles_InP_small

# Enthalpy contribution from oleate protonation (per mole of InP)
dH_from_oleate_small = molar_ratio_small * dH_prot_oleate

# --- Final Analysis ---
# Calculate the change in measured dissolution enthalpy
observed_enthalpy_change = dH_diss_small_measured - dH_diss_large_measured

# Calculate the change in enthalpy attributable to oleate protonation
protonation_enthalpy_change = dH_from_oleate_small - dH_from_oleate_large

# Display the final equation and results
# Using sys.stdout.write to prevent the final answer tag from being captured by the code interpreter.
sys.stdout.write("To evaluate hypothesis A, we compare the observed enthalpy change to the change from oleate protonation.\n\n")
sys.stdout.write("Observed Enthalpy Change = ΔH_small - ΔH_large\n")
sys.stdout.write(f"Observed Enthalpy Change = {dH_diss_small_measured} kJ/mol - {dH_diss_large_measured} kJ/mol = {observed_enthalpy_change:.2f} kJ/mol\n\n")
sys.stdout.write("Calculated Change from Protonation = (ΔH from Oleate)_small - (ΔH from Oleate)_large\n")
sys.stdout.write(f"Calculated Change from Protonation = {dH_from_oleate_small:.2f} kJ/mol - {dH_from_oleate_large:.2f} kJ/mol = {protonation_enthalpy_change:.2f} kJ/mol\n\n")
sys.stdout.write(f"Conclusion: The enthalpy change from oleate protonation ({protonation_enthalpy_change:.2f} kJ/mol) only accounts for a very small fraction of the total observed change ({observed_enthalpy_change:.2f} kJ/mol). Therefore, hypothesis A is not the primary explanation.\n")
sys.stdout.write("Among the choices, option D provides the most logical alternative explanation for the significant remaining endothermic contribution.\n")
