# --- Constants and Given Data ---
# Molecular weights (g/mol)
MW_INP = 146
MW_OLEATE = 281

# Enthalpies (kJ/mol)
DH_DISSOLUTION_BULK_INP = -86
DH_PROTONATION_OLEATE = 7

# Data for the largest quantum dot (LQD)
MASS_FRAC_OLEATE_LQD = 0.46  # 46%
DH_DISSOLUTION_LQD = 70      # kJ/mol of InP

# Data for the smallest quantum dot (SQD)
MASS_FRAC_OLEATE_SQD = 0.52  # 52%
DH_DISSOLUTION_SQD = 120     # kJ/mol of InP

# --- Calculation for the Largest Quantum Dot (LQD) ---
# Assume a 100g sample
mass_inp_lqd = 100 * (1 - MASS_FRAC_OLEATE_LQD)
mass_oleate_lqd = 100 * MASS_FRAC_OLEATE_LQD

# Convert mass to moles
moles_inp_lqd = mass_inp_lqd / MW_INP
moles_oleate_lqd = mass_oleate_lqd / MW_OLEATE

# Calculate mole ratio of oleate to InP
ratio_lqd = moles_oleate_lqd / moles_inp_lqd

# Calculate enthalpy contribution from oleate protonation per mole of InP
contrib_oleate_lqd = ratio_lqd * DH_PROTONATION_OLEATE

# --- Calculation for the Smallest Quantum Dot (SQD) ---
# Assume a 100g sample
mass_inp_sqd = 100 * (1 - MASS_FRAC_OLEATE_SQD)
mass_oleate_sqd = 100 * MASS_FRAC_OLEATE_SQD

# Convert mass to moles
moles_inp_sqd = mass_inp_sqd / MW_INP
moles_oleate_sqd = mass_oleate_sqd / MW_OLEATE

# Calculate mole ratio of oleate to InP
ratio_sqd = moles_oleate_sqd / moles_inp_sqd

# Calculate enthalpy contribution from oleate protonation per mole of InP
contrib_oleate_sqd = ratio_sqd * DH_PROTONATION_OLEATE

# --- Analysis of Results ---
# Total observed change in enthalpy
total_dh_change = DH_DISSOLUTION_SQD - DH_DISSOLUTION_LQD

# Change in enthalpy due to oleate protonation
dh_change_from_oleate = contrib_oleate_sqd - contrib_oleate_lqd

# --- Print the step-by-step analysis ---
print("--- Analysis of Oleate Protonation Contribution ---")
print("\nStep 1: Calculate moles of oleate per mole of InP for the Largest QD")
print(f"For every {moles_inp_lqd:.3f} moles of InP, there are {moles_oleate_lqd:.3f} moles of oleate.")
print(f"Mole Ratio (oleate/InP) = {moles_oleate_lqd:.3f} / {moles_inp_lqd:.3f} = {ratio_lqd:.3f}")
print("\nStep 2: Calculate enthalpy contribution from oleate protonation for the Largest QD")
print(f"Contribution = Mole Ratio * Enthalpy of Protonation")
print(f"Contribution = {ratio_lqd:.3f} * {DH_PROTONATION_OLEATE} kJ/mol = {contrib_oleate_lqd:.2f} kJ/mol of InP")

print("\n-----------------------------------------------------")

print("\nStep 3: Calculate moles of oleate per mole of InP for the Smallest QD")
print(f"For every {moles_inp_sqd:.3f} moles of InP, there are {moles_oleate_sqd:.3f} moles of oleate.")
print(f"Mole Ratio (oleate/InP) = {moles_oleate_sqd:.3f} / {moles_inp_sqd:.3f} = {ratio_sqd:.3f}")
print("\nStep 4: Calculate enthalpy contribution from oleate protonation for the Smallest QD")
print(f"Contribution = Mole Ratio * Enthalpy of Protonation")
print(f"Contribution = {ratio_sqd:.3f} * {DH_PROTONATION_OLEATE} kJ/mol = {contrib_oleate_sqd:.2f} kJ/mol of InP")

print("\n-----------------------------------------------------")
print("\n--- Final Comparison ---")
print(f"The total observed increase in dissolution enthalpy is:")
print(f"{DH_DISSOLUTION_SQD} kJ/mol - {DH_DISSOLUTION_LQD} kJ/mol = {total_dh_change} kJ/mol")
print("\nThe increase in enthalpy that can be explained by oleate protonation is:")
print(f"{contrib_oleate_sqd:.2f} kJ/mol - {contrib_oleate_lqd:.2f} kJ/mol = {dh_change_from_oleate:.2f} kJ/mol")
print("\nConclusion: The change due to oleate protonation (~0.8 kJ/mol) is a very small fraction of the total observed change (50 kJ/mol). Therefore, the increase in oleate amount and its protonation is not the primary explanation for the large endothermic shift.")
