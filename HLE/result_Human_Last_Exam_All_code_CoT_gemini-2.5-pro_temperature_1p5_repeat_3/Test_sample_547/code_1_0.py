import math

# --- Given Constants ---
MW_InP = 146  # g/mol
MW_oleate = 281  # g/mol
dH_diss_bulk_InP = -86  # kJ/mol
dH_protonation_oleate = 7  # kJ/mol

# --- Data for Large Quantum Dot (LQD) ---
mass_frac_oleate_LQD = 0.46
dH_diss_total_LQD = 70  # kJ/mol of InP

# --- Data for Small Quantum Dot (SQD) ---
mass_frac_oleate_SQD = 0.52
dH_diss_total_SQD = 120  # kJ/mol of InP

# --- Calculations for a hypothetical 100g sample ---

# For Large QD
mass_InP_LQD = 100 * (1 - mass_frac_oleate_LQD)
mass_oleate_LQD = 100 * mass_frac_oleate_LQD
moles_InP_LQD = mass_InP_LQD / MW_InP
moles_oleate_LQD = mass_oleate_LQD / MW_oleate
ratio_LQD = moles_oleate_LQD / moles_InP_LQD
dH_protonation_contribution_LQD = ratio_LQD * dH_protonation_oleate

# For Small QD
mass_InP_SQD = 100 * (1 - mass_frac_oleate_SQD)
mass_oleate_SQD = 100 * mass_frac_oleate_SQD
moles_InP_SQD = mass_InP_SQD / MW_InP
moles_oleate_SQD = mass_oleate_SQD / MW_oleate
ratio_SQD = moles_oleate_SQD / moles_InP_SQD
dH_protonation_contribution_SQD = ratio_SQD * dH_protonation_oleate

# --- Analysis ---
total_dH_change = dH_diss_total_SQD - dH_diss_total_LQD
protonation_dH_change = dH_protonation_contribution_SQD - dH_protonation_contribution_LQD
percentage_explained_by_protonation = (protonation_dH_change / total_dH_change) * 100

# Enthalpy not accounted for by bulk dissolution or protonation
# This is the "excess surface enthalpy" which includes InP surface energy and ligand shell effects
excess_enthalpy_LQD = dH_diss_total_LQD - dH_diss_bulk_InP - dH_protonation_contribution_LQD
excess_enthalpy_SQD = dH_diss_total_SQD - dH_diss_bulk_InP - dH_protonation_contribution_SQD

# --- Output the results ---
print("--- Analysis of Enthalpy Contributions ---")
print("\n# Large Quantum Dot:")
print(f"Mole ratio (oleate/InP): {ratio_LQD:.4f}")
print(f"Enthalpy from oleate protonation (per mol InP): {dH_protonation_contribution_LQD:.2f} kJ/mol")
print(f"Total dissolution enthalpy: {dH_diss_total_LQD} kJ/mol")
print(f"Contribution from bulk InP dissolution: {dH_diss_bulk_InP} kJ/mol")
print(f"Combined excess surface enthalpy: {excess_enthalpy_LQD:.2f} kJ/mol = {dH_diss_total_LQD} - ({dH_diss_bulk_InP}) - {dH_protonation_contribution_LQD:.2f}")

print("\n# Small Quantum Dot:")
print(f"Mole ratio (oleate/InP): {ratio_SQD:.4f}")
print(f"Enthalpy from oleate protonation (per mol InP): {dH_protonation_contribution_SQD:.2f} kJ/mol")
print(f"Total dissolution enthalpy: {dH_diss_total_SQD} kJ/mol")
print(f"Contribution from bulk InP dissolution: {dH_diss_bulk_InP} kJ/mol")
print(f"Combined excess surface enthalpy: {excess_enthalpy_SQD:.2f} kJ/mol = {dH_diss_total_SQD} - ({dH_diss_bulk_InP}) - {dH_protonation_contribution_SQD:.2f}")


print("\n# Comparison and Conclusion:")
print(f"Total observed change in enthalpy (Small QD - Large QD): {total_dH_change:.2f} kJ/mol = {dH_diss_total_SQD} - {dH_diss_total_LQD}")
print(f"Change due to oleate protonation: {protonation_dH_change:.2f} kJ/mol = {dH_protonation_contribution_SQD:.2f} - {dH_protonation_contribution_LQD:.2f}")
print(f"Percentage of total change explained by protonation alone: {percentage_explained_by_protonation:.2f}%")
print("\nThis shows that the endothermic protonation of the increased amount of oleate (Choice A) only accounts for a very small fraction of the observed increase in enthalpy.")
print("The majority of the change must come from other surface-related effects that are more pronounced in smaller dots, such as the energy required to disrupt the ordered ligand shell (Choice D) or the surface energy of the InP core itself.")
print("Choice C is incorrect because surface energy is positive, not negative. Choice D provides a physically sound explanation for a significant endothermic contribution that scales with the amount of ligand.")
