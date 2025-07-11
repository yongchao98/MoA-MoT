import pandas as pd

# --- Given Constants ---
MW_InP = 146  # g/mol
MW_oleate = 281  # g/mol
dH_protonation_oleate = 7  # kJ/mol of oleate

# --- Data for the Quantum Dots ---
# Data is structured for easy calculation
qd_data = {
    'Largest QD': {
        'mass_fraction_oleate': 0.46,
        'dH_dissolution_total': 70  # kJ/mol of InP
    },
    'Smallest QD': {
        'mass_fraction_oleate': 0.52,
        'dH_dissolution_total': 120  # kJ/mol of InP
    }
}

# --- Calculation ---
results = {}
basis_mass = 100  # Assume a 100g sample for calculation

for name, data in qd_data.items():
    # Mass fractions
    mass_frac_oleate = data['mass_fraction_oleate']
    mass_frac_InP = 1 - mass_frac_oleate
    
    # Mass in the 100g sample
    mass_oleate = basis_mass * mass_frac_oleate
    mass_InP = basis_mass * mass_frac_InP
    
    # Moles in the 100g sample
    moles_oleate = mass_oleate / MW_oleate
    moles_InP = mass_InP / MW_InP
    
    # Molar ratio of oleate to InP
    ratio_oleate_to_InP = moles_oleate / moles_InP
    
    # Enthalpy contribution from oleate protonation, per mole of InP
    dH_from_protonation = ratio_oleate_to_InP * dH_protonation_oleate
    
    # Store results
    results[name] = {
        'moles_oleate_per_mole_InP': ratio_oleate_to_InP,
        'dH_from_protonation (kJ/mol InP)': dH_from_protonation
    }

# --- Analysis of the Change ---
# Observed total change in dissolution enthalpy
total_dH_change = qd_data['Smallest QD']['dH_dissolution_total'] - qd_data['Largest QD']['dH_dissolution_total']

# Change in enthalpy attributed to oleate protonation
protonation_dH_change = results['Smallest QD']['dH_from_protonation (kJ/mol InP)'] - results['Largest QD']['dH_from_protonation (kJ/mol InP)']

# --- Print Results ---
print("--- Calculation Results ---")
df = pd.DataFrame(results).T
print(df.to_string())
print("\n--- Analysis of Enthalpy Change ---")
print(f"The observed total change in dissolution enthalpy is: {qd_data['Largest QD']['dH_dissolution_total']} -> {qd_data['Smallest QD']['dH_dissolution_total']} kJ/mol InP.")
print(f"Total increase in endothermicity = {total_dH_change:.1f} kJ/mol InP.")
print(f"The calculated change in enthalpy from oleate protonation is: {results['Largest QD']['dH_from_protonation (kJ/mol InP)']:0.2f} -> {results['Smallest QD']['dH_from_protonation (kJ/mol InP)']:0.2f} kJ/mol InP.")
print(f"Increase in endothermicity from oleate protonation alone = {protonation_dH_change:.2f} kJ/mol InP.")
print(f"\nThe oleate protonation only accounts for {protonation_dH_change:.2f} / {total_dH_change:.1f} = {protonation_dH_change/total_dH_change:.1%} of the total observed change.")

print("\n--- Conclusion ---")
print("The calculation shows that while the endothermic protonation of oleate contributes to the overall enthalpy, its increase from the largest to the smallest quantum dot is only ~0.82 kJ/mol.")
print("This is a very small fraction of the total observed increase of 50 kJ/mol.")
print("Therefore, answer A is not the most logical explanation for the large observed change.")
print("Answer D provides a more plausible explanation: disrupting the dense ligand shell requires a significant amount of energy, and this effect becomes larger for smaller dots which have a greater proportion of ligands. The energy required to break these strong inter-ligand (van der Waals) interactions is substantial and scales correctly with the observations.")
