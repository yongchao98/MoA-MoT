# Given constants
mw_inp = 146  # g/mol
mw_oleate = 281  # g/mol
h_protonation_oleate = 7  # kJ/mol of oleate

# Data for the largest quantum dot (QD)
h_diss_total_large = 70  # kJ/mol of InP
mass_frac_oleate_large = 0.46

# Data for the smallest quantum dot (QD)
h_diss_total_small = 120  # kJ/mol of InP
mass_frac_oleate_small = 0.52

print("Analyzing the enthalpy contributions for the dissolution of quantum dots.")
print("The total enthalpy equation is: ΔH_total = ΔH_oleate + ΔH_other\n")

# --- Calculations for the LARGEST quantum dot ---

# Assume 100g of the QD material
mass_oleate_large = 100 * mass_frac_oleate_large
mass_inp_large = 100 * (1 - mass_frac_oleate_large)

# Moles of each component
moles_oleate_large = mass_oleate_large / mw_oleate
moles_inp_large = mass_inp_large / mw_inp

# Molar ratio of oleate to InP
molar_ratio_large = moles_oleate_large / moles_inp_large

# Enthalpy contribution from oleate protonation (per mole of InP)
h_oleate_large = molar_ratio_large * h_protonation_oleate

# The rest of the enthalpy contribution
h_other_large = h_diss_total_large - h_oleate_large

print("--- For the LARGEST Quantum Dot ---")
print(f"Molar ratio of Oleate to InP: {molar_ratio_large:.3f}")
print(f"Enthalpy from Oleate Protonation (ΔH_oleate) = {h_oleate_large:.2f} kJ/mol")
print(f"Enthalpy from other sources (ΔH_other) = {h_other_large:.2f} kJ/mol")
print("Full equation:")
print(f"{h_diss_total_large:.2f} (Total) = {h_oleate_large:.2f} (Oleate) + {h_other_large:.2f} (Other)\n")


# --- Calculations for the SMALLEST quantum dot ---

# Assume 100g of the QD material
mass_oleate_small = 100 * mass_frac_oleate_small
mass_inp_small = 100 * (1 - mass_frac_oleate_small)

# Moles of each component
moles_oleate_small = mass_oleate_small / mw_oleate
moles_inp_small = mass_inp_small / mw_inp

# Molar ratio of oleate to InP
molar_ratio_small = moles_oleate_small / moles_inp_small

# Enthalpy contribution from oleate protonation (per mole of InP)
h_oleate_small = molar_ratio_small * h_protonation_oleate

# The rest of the enthalpy contribution
h_other_small = h_diss_total_small - h_oleate_small

print("--- For the SMALLEST Quantum Dot ---")
print(f"Molar ratio of Oleate to InP: {molar_ratio_small:.3f}")
print(f"Enthalpy from Oleate Protonation (ΔH_oleate) = {h_oleate_small:.2f} kJ/mol")
print(f"Enthalpy from other sources (ΔH_other) = {h_other_small:.2f} kJ/mol")
print("Full equation:")
print(f"{h_diss_total_small:.2f} (Total) = {h_oleate_small:.2f} (Oleate) + {h_other_small:.2f} (Other)\n")


# --- Comparison and Conclusion ---
total_h_change = h_diss_total_small - h_diss_total_large
oleate_h_change = h_oleate_small - h_oleate_large

print("--- Analysis of Change ---")
print(f"The total measured enthalpy becomes more endothermic by: {total_h_change:.2f} kJ/mol")
print(f"The change due to oleate protonation is only: {oleate_h_change:.2f} kJ/mol")
print("\nConclusion: The increased amount of oleate (Answer A) only accounts for a tiny fraction of the total change in enthalpy.")
print("The vast majority of the change comes from other size-dependent effects, such as the energy required to disrupt the packed ligand shell.")
