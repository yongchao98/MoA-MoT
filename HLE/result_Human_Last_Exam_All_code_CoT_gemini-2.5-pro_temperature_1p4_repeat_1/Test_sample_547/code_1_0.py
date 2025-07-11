# Given data
mw_inp = 146  # g/mol
mw_oleate = 281  # g/mol
h_protonation_oleate = 7  # kJ/mol

# Largest quantum dot data
mass_frac_oleate_large = 0.46
h_dissolution_large = 70  # kJ/mol

# Smallest quantum dot data
mass_frac_oleate_small = 0.52
h_dissolution_small = 120  # kJ/mol

# --- Calculations for the largest QD ---
# Consider a 100g sample
mass_oleate_large = 100 * mass_frac_oleate_large
mass_inp_large = 100 * (1 - mass_frac_oleate_large)

# Moles in the 100g sample
moles_oleate_large = mass_oleate_large / mw_oleate
moles_inp_large = mass_inp_large / mw_inp

# Molar ratio of oleate to InP
ratio_large = moles_oleate_large / moles_inp_large

# Enthalpy contribution from oleate protonation per mole of InP
h_oleate_contrib_large = ratio_large * h_protonation_oleate

# --- Calculations for the smallest QD ---
# Consider a 100g sample
mass_oleate_small = 100 * mass_frac_oleate_small
mass_inp_small = 100 * (1 - mass_frac_oleate_small)

# Moles in the 100g sample
moles_oleate_small = mass_oleate_small / mw_oleate
moles_inp_small = mass_inp_small / mw_inp

# Molar ratio of oleate to InP
ratio_small = moles_oleate_small / moles_inp_small

# Enthalpy contribution from oleate protonation per mole of InP
h_oleate_contrib_small = ratio_small * h_protonation_oleate

# --- Compare calculated and observed changes ---
# Calculated change in enthalpy due to oleate protonation
delta_h_oleate_calc = h_oleate_contrib_small - h_oleate_contrib_large

# Observed total change in dissolution enthalpy
delta_h_observed = h_dissolution_small - h_dissolution_large

# --- Print results ---
print("Analysis of Enthalpy Contributions:")
print("-" * 40)
print("For the Largest Quantum Dot:")
print(f"The molar ratio of oleate to InP is: {ratio_large:.4f}")
print(f"Enthalpy contribution from oleate protonation: {ratio_large:.4f} * {h_protonation_oleate} kJ/mol = {h_oleate_contrib_large:.2f} kJ/mol of InP")
print()
print("For the Smallest Quantum Dot:")
print(f"The molar ratio of oleate to InP is: {ratio_small:.4f}")
print(f"Enthalpy contribution from oleate protonation: {ratio_small:.4f} * {h_protonation_oleate} kJ/mol = {h_oleate_contrib_small:.2f} kJ/mol of InP")
print()
print("Comparison:")
print(f"Observed total increase in endothermic enthalpy: {h_dissolution_small} - {h_dissolution_large} = {delta_h_observed:.2f} kJ/mol")
print(f"Calculated increase from oleate protonation alone: {h_oleate_contrib_small:.2f} - {h_oleate_contrib_large:.2f} = {delta_h_oleate_calc:.2f} kJ/mol")
print("-" * 40)
print("\nConclusion:")
print("The increased amount of oleate in smaller quantum dots only accounts for a very small fraction")
print(f"({delta_h_oleate_calc:.2f} kJ/mol) of the total observed enthalpy change ({delta_h_observed:.2f} kJ/mol).")
print("Therefore, another more significant energy contribution must be responsible.")