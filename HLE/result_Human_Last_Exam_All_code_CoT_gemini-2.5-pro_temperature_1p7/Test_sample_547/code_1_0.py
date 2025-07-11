# --- Input Data ---
# Molecular Weights (g/mol)
MW_InP = 146
MW_oleate = 281

# Enthalpies (kJ/mol)
dH_protonation_oleate = 7
dH_diss_large_QD = 70
dH_diss_small_QD = 120

# Mass Fractions
mass_fraction_oleate_large = 0.46
mass_fraction_oleate_small = 0.52

# --- Step 1: Calculate molar ratio for the large quantum dot ---
# For a conceptual 100g sample of large QDs:
# Mass of InP = 100g * (1 - 0.46) = 54g
# Mass of Oleate = 100g * 0.46 = 46g
# Moles of InP = 54g / 146 g/mol
# Moles of Oleate = 46g / 281 g/mol
moles_InP_large_per_100g = (1 - mass_fraction_oleate_large) * 100 / MW_InP
moles_oleate_large_per_100g = mass_fraction_oleate_large * 100 / MW_oleate
molar_ratio_large = moles_oleate_large_per_100g / moles_InP_large_per_100g

# --- Step 2: Calculate protonation enthalpy contribution for the large quantum dot ---
dH_protonation_contribution_large = molar_ratio_large * dH_protonation_oleate

# --- Step 1 (cont.): Calculate molar ratio for the small quantum dot ---
# For a conceptual 100g sample of small QDs:
# Mass of InP = 100g * (1 - 0.52) = 48g
# Mass of Oleate = 100g * 0.52 = 52g
# Moles of InP = 48g / 146 g/mol
# Moles of Oleate = 52g / 281 g/mol
moles_InP_small_per_100g = (1 - mass_fraction_oleate_small) * 100 / MW_InP
moles_oleate_small_per_100g = mass_fraction_oleate_small * 100 / MW_oleate
molar_ratio_small = moles_oleate_small_per_100g / moles_InP_small_per_100g

# --- Step 2 (cont.): Calculate protonation enthalpy contribution for the small quantum dot ---
dH_protonation_contribution_small = molar_ratio_small * dH_protonation_oleate

# --- Step 3 & 4: Compare the changes and analyze ---
# Total observed change in enthalpy
total_dH_change = dH_diss_small_QD - dH_diss_large_QD
# Change in enthalpy due to oleate protonation
protonation_dH_change = dH_protonation_contribution_small - dH_protonation_contribution_large

print("Analysis of Enthalpy Change in Quantum Dot Dissolution")
print("="*55)
print("The overall observed change in dissolution enthalpy is:")
print(f"{dH_diss_small_QD} kJ/mol (small QD) - {dH_diss_large_QD} kJ/mol (large QD) = {total_dH_change:.2f} kJ/mol\n")

print("Let's calculate the portion of this change caused by oleate protonation (Choice A).")
print(f"For the large QD, the protonation of oleate contributes {molar_ratio_large:.2f} * {dH_protonation_oleate} = {dH_protonation_contribution_large:.2f} kJ per mole of InP.")
print(f"For the small QD, the protonation of oleate contributes {molar_ratio_small:.2f} * {dH_protonation_oleate} = {dH_protonation_contribution_small:.2f} kJ per mole of InP.\n")

print("The difference in enthalpy due to increased oleate is therefore:")
print(f"{dH_protonation_contribution_small:.2f} kJ/mol - {dH_protonation_contribution_large:.2f} kJ/mol = {protonation_dH_change:.2f} kJ/mol\n")

print("Conclusion:")
print(f"The increased amount of oleate only accounts for {protonation_dH_change:.2f} kJ/mol of the {total_dH_change:.2f} kJ/mol total endothermic shift.")
print("This contribution is too small to be the main explanation.")
print("Therefore, Choice A is not the most logical explanation.")
print("Choice D proposes that disrupting the packed ligand shell is an endothermic process that becomes larger for smaller dots with more ligands, which better explains the large magnitude of the observed enthalpy change.")