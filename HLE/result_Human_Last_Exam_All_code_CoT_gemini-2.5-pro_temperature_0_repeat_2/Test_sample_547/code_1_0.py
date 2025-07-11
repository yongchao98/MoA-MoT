# Define constants from the problem description
MW_InP = 146  # g/mol
MW_oleate = 281  # g/mol

# Data for the largest quantum dot
w_oleate_large = 0.46  # 46%
H_diss_large = 70  # kJ/mol of InP

# Data for the smallest quantum dot
w_oleate_small = 0.52  # 52%
H_diss_small = 120  # kJ/mol of InP

# Enthalpy of protonation for oleate
H_protonation_oleate = 7  # kJ/mol

# --- Calculations for the largest quantum dot ---
# Assume a 100g sample
mass_oleate_large = 100 * w_oleate_large
mass_InP_large = 100 * (1 - w_oleate_large)

# Calculate moles
moles_oleate_large = mass_oleate_large / MW_oleate
moles_InP_large = mass_InP_large / MW_InP

# Calculate mole ratio of oleate to InP
mole_ratio_large = moles_oleate_large / moles_InP_large

# Calculate the enthalpy contribution from oleate protonation per mole of InP
H_contrib_oleate_large = mole_ratio_large * H_protonation_oleate

# --- Calculations for the smallest quantum dot ---
# Assume a 100g sample
mass_oleate_small = 100 * w_oleate_small
mass_InP_small = 100 * (1 - w_oleate_small)

# Calculate moles
moles_oleate_small = mass_oleate_small / MW_oleate
moles_InP_small = mass_InP_small / MW_InP

# Calculate mole ratio of oleate to InP
mole_ratio_small = moles_oleate_small / moles_InP_small

# Calculate the enthalpy contribution from oleate protonation per mole of InP
H_contrib_oleate_small = mole_ratio_small * H_protonation_oleate

# --- Analysis ---
# Calculate the change in enthalpy due to oleate protonation
delta_H_oleate = H_contrib_oleate_small - H_contrib_oleate_large

# Calculate the total observed change in dissolution enthalpy
delta_H_observed = H_diss_small - H_diss_large

# --- Output the results ---
print("Analysis of Enthalpy Contributions:")
print("-" * 40)
print("For the LARGEST quantum dot:")
print(f"The mole ratio of oleate to InP is: {mole_ratio_large:.3f}")
print(f"The enthalpy contribution from oleate protonation is {mole_ratio_large:.3f} * {H_protonation_oleate} kJ/mol = {H_contrib_oleate_large:.2f} kJ/mol of InP.")
print("-" * 40)
print("For the SMALLEST quantum dot:")
print(f"The mole ratio of oleate to InP is: {mole_ratio_small:.3f}")
print(f"The enthalpy contribution from oleate protonation is {mole_ratio_small:.3f} * {H_protonation_oleate} kJ/mol = {H_contrib_oleate_small:.2f} kJ/mol of InP.")
print("-" * 40)
print("Comparison:")
print(f"The total observed increase in endothermic enthalpy is {H_diss_small} - {H_diss_large} = {delta_H_observed:.2f} kJ/mol of InP.")
print(f"The calculated increase in enthalpy from oleate protonation is {H_contrib_oleate_small:.2f} - {H_contrib_oleate_large:.2f} = {delta_H_oleate:.2f} kJ/mol of InP.")
print("\nConclusion:")
print(f"The increased amount of oleate only accounts for {delta_H_oleate:.2f} kJ/mol of the total {delta_H_observed:.2f} kJ/mol change.")
print("This is a very small fraction of the observed effect. Therefore, another, more significant endothermic process that scales with surface area must be responsible.")
print("This supports choice D, which posits that disrupting the tightly packed ligand shell is a significant endothermic process.")
