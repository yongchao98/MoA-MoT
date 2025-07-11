import math

# --- Given constants and data ---
MW_InP = 146  # g/mol
MW_oleate = 281  # g/mol
dH_protonation_oleate = 7  # kJ/mol

# Data for the largest quantum dot
mass_frac_oleate_large = 0.46  # 46%
dH_diss_large = 70  # kJ/mol of InP

# Data for the smallest quantum dot
mass_frac_oleate_small = 0.52  # 52%
dH_diss_small = 120  # kJ/mol of InP

def calculate_protonation_contribution(mass_frac_oleate):
    """
    Calculates the enthalpy contribution from oleate protonation per mole of InP.
    We assume a 100g sample for calculation.
    """
    mass_oleate = 100 * mass_frac_oleate
    mass_InP = 100 * (1 - mass_frac_oleate)
    
    moles_oleate = mass_oleate / MW_oleate
    moles_InP = mass_InP / MW_InP
    
    # Molar ratio of oleate to InP
    moles_oleate_per_mol_InP = moles_oleate / moles_InP
    
    # Enthalpy contribution from oleate protonation per mole of InP
    enthalpy_contribution = moles_oleate_per_mol_InP * dH_protonation_oleate
    
    return enthalpy_contribution

# --- Calculations ---
# Contribution for the largest QD
protonation_H_large = calculate_protonation_contribution(mass_frac_oleate_large)

# Contribution for the smallest QD
protonation_H_small = calculate_protonation_contribution(mass_frac_oleate_small)

# Change in enthalpy due to protonation
protonation_H_change = protonation_H_small - protonation_H_large

# Observed change in total enthalpy
observed_H_change = dH_diss_small - dH_diss_large


# --- Output the results ---
print("--- Analysis of Enthalpy Contributions ---")
print(f"For the largest QD, the enthalpy from oleate protonation is: {protonation_H_large:.2f} kJ/mol of InP")
print(f"For the smallest QD, the enthalpy from oleate protonation is: {protonation_H_small:.2f} kJ/mol of InP")
print("\n--- Comparing the change in enthalpy ---")
print(f"The observed increase in dissolution enthalpy is: {dH_diss_small} kJ/mol - {dH_diss_large} kJ/mol = {observed_H_change:.2f} kJ/mol")
print(f"The calculated increase from oleate protonation alone is: {protonation_H_small:.2f} kJ/mol - {protonation_H_large:.2f} kJ/mol = {protonation_H_change:.2f} kJ/mol")

print("\n--- Conclusion ---")
print("The contribution from oleate protonation accounts for only a very small fraction of the total observed change in enthalpy.")
print("Therefore, another, more significant endothermic process that scales with the surface area must be responsible.")
print("This points to the energy required to break apart the stabilized ligand shell (Choice D) as the most logical explanation.")
