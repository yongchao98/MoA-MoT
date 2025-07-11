# Constants and given data
MW_InP = 146  # g/mol
MW_oleate = 281  # g/mol
dH_diss_bulk_InP = -86  # kJ/mol of InP
dH_protonation_oleate = 7  # kJ/mol of oleate

# Data for the largest quantum dot
mass_frac_oleate_large = 0.46
dH_diss_large_QD = 70  # kJ/mol of InP

# Data for the smallest quantum dot
mass_frac_oleate_small = 0.52
dH_diss_small_QD = 120  # kJ/mol of InP

print("Step 1: Calculate the enthalpy contribution from oleate protonation.")

# --- Largest Quantum Dot Calculation ---
# Based on a 100g sample of QDs
mass_InP_large = 100 * (1 - mass_frac_oleate_large)
mass_oleate_large = 100 * mass_frac_oleate_large
moles_InP_large = mass_InP_large / MW_InP
moles_oleate_large = mass_oleate_large / MW_oleate
molar_ratio_large = moles_oleate_large / moles_InP_large
dH_oleate_contrib_large = molar_ratio_large * dH_protonation_oleate

print("\n--- For the Largest Quantum Dot ---")
print(f"The molar ratio of oleate to InP is {molar_ratio_large:.4f}.")
print("The enthalpy from oleate protonation is calculated as:")
print(f"Contribution = (molar ratio) * (protonation enthalpy of oleate)")
print(f"Contribution = {molar_ratio_large:.4f} mol_oleate/mol_InP * {dH_protonation_oleate} kJ/mol_oleate = {dH_oleate_contrib_large:.2f} kJ/mol_InP")


# --- Smallest Quantum Dot Calculation ---
# Based on a 100g sample of QDs
mass_InP_small = 100 * (1 - mass_frac_oleate_small)
mass_oleate_small = 100 * mass_frac_oleate_small
moles_InP_small = mass_InP_small / MW_InP
moles_oleate_small = mass_oleate_small / MW_oleate
molar_ratio_small = moles_oleate_small / moles_InP_small
dH_oleate_contrib_small = molar_ratio_small * dH_protonation_oleate

print("\n--- For the Smallest Quantum Dot ---")
print(f"The molar ratio of oleate to InP is {molar_ratio_small:.4f}.")
print("The enthalpy from oleate protonation is calculated as:")
print(f"Contribution = {molar_ratio_small:.4f} mol_oleate/mol_InP * {dH_protonation_oleate} kJ/mol_oleate = {dH_oleate_contrib_small:.2f} kJ/mol_InP")

print("\nStep 2: Determine if oleate protonation explains the observed enthalpy trend.")
total_dH_change_observed = dH_diss_small_QD - dH_diss_large_QD
dH_change_from_oleate = dH_oleate_contrib_small - dH_oleate_contrib_large
print(f"The observed enthalpy increase from large to small QDs is: {dH_diss_small_QD} - {dH_diss_large_QD} = {total_dH_change_observed} kJ/mol.")
print(f"The increase that can be explained by oleate protonation is only: {dH_oleate_contrib_small:.2f} - {dH_oleate_contrib_large:.2f} = {dH_change_from_oleate:.2f} kJ/mol.")
print("This contribution is far too small. Therefore, option A is incorrect.")


print("\nStep 3: Calculate the large, unexplained enthalpy term.")
print("We can model the total enthalpy as: ΔH_total = ΔH_bulk_InP + ΔH_oleate + ΔH_other")
print("Therefore: ΔH_other = ΔH_total - ΔH_bulk_InP - ΔH_oleate")

# Calculation for the 'other' enthalpy term for the largest QD
dH_other_large = dH_diss_large_QD - dH_diss_bulk_InP - dH_oleate_contrib_large
print(f"\nFor the largest QD, the unexplained enthalpy (ΔH_other) is:")
print(f"ΔH_other = {dH_diss_large_QD} - ({dH_diss_bulk_InP}) - {dH_oleate_contrib_large:.2f} = {dH_other_large:.2f} kJ/mol")

# Calculation for the 'other' enthalpy term for the smallest QD
dH_other_small = dH_diss_small_QD - dH_diss_bulk_InP - dH_oleate_contrib_small
print(f"\nFor the smallest QD, the unexplained enthalpy (ΔH_other) is:")
print(f"ΔH_other = {dH_diss_small_QD} - ({dH_diss_bulk_InP}) - {dH_oleate_contrib_small:.2f} = {dH_other_small:.2f} kJ/mol")

print("\nFinal Conclusion:")
print("The calculation reveals a very large endothermic term (152.9 to 202.1 kJ/mol) that is not explained by the dissolution of the InP core or ligand protonation.")
print("This term increases as the quantum dot size decreases, which correlates with an increasing proportion of ligands.")
print("Option D, which attributes this energy cost to disrupting the tightly packed shell of interacting organic ligands, provides the most logical explanation for this large, size-dependent endothermic effect.")
