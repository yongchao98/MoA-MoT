import sys
# Redirect stdout to a string buffer to capture all prints
from io import StringIO
old_stdout = sys.stdout
sys.stdout = captured_output = StringIO()


# --- Given Constants ---
MW_InP = 146  # g/mol
MW_oleate = 281  # g/mol
dH_bulk_InP = -86  # kJ/mol
dH_oleate_protonation = 7  # kJ/mol

# --- Data for the Largest Quantum Dot ---
mass_fraction_oleate_large = 0.46  # 46%
dH_measured_large = 70  # kJ/mol

# --- Data for the Smallest Quantum Dot ---
mass_fraction_oleate_small = 0.52  # 52%
dH_measured_small = 120  # kJ/mol

def analyze_quantum_dot(name, mass_fraction_oleate, dH_measured):
    """
    Calculates the theoretical and missing enthalpy for a given quantum dot.
    """
    print(f"--- Analysis for the {name} Quantum Dot ---")
    
    # Step 1: Calculate molar ratio of oleate to InP
    # Assume a 100g sample for calculation
    mass_oleate = mass_fraction_oleate * 100
    mass_inp = (1 - mass_fraction_oleate) * 100
    
    moles_oleate = mass_oleate / MW_oleate
    moles_inp = mass_inp / MW_InP
    
    molar_ratio_oleate_to_inp = moles_oleate / moles_inp
    
    print(f"Step 1: Calculate the molar ratio of oleate to InP.")
    print(f"For a mass fraction of {mass_fraction_oleate*100:.0f}% oleate, the molar ratio of oleate to InP is {molar_ratio_oleate_to_inp:.3f}.")
    print("-" * 20)
    
    # Step 2: Calculate enthalpy contribution from oleate protonation
    dH_oleate_contribution = molar_ratio_oleate_to_inp * dH_oleate_protonation
    
    print(f"Step 2: Calculate the enthalpy from oleate protonation per mole of InP.")
    print(f"Enthalpy = (Molar Ratio) * (Oleate Protonation Enthalpy)")
    print(f"Enthalpy = {molar_ratio_oleate_to_inp:.3f} * {dH_oleate_protonation} kJ/mol = {dH_oleate_contribution:.2f} kJ/mol of InP.")
    print("-" * 20)

    # Step 3: Calculate the total expected enthalpy based on bulk values
    dH_expected = dH_bulk_InP + dH_oleate_contribution
    
    print(f"Step 3: Calculate the total expected enthalpy from core dissolution and ligand protonation.")
    print(f"Expected Enthalpy = (Bulk InP Enthalpy) + (Oleate Contribution)")
    print(f"Expected Enthalpy = {dH_bulk_InP} kJ/mol + {dH_oleate_contribution:.2f} kJ/mol = {dH_expected:.2f} kJ/mol.")
    print("-" * 20)

    # Step 4: Compare with the measured value and find the discrepancy
    missing_enthalpy = dH_measured - dH_expected
    
    print(f"Step 4: Compare the expected enthalpy with the measured value.")
    print(f"The measured enthalpy is {dH_measured} kJ/mol.")
    print(f"The calculated expected enthalpy is only {dH_expected:.2f} kJ/mol.")
    print(f"Discrepancy (Missing Enthalpy) = {dH_measured} - ({dH_expected:.2f}) = {missing_enthalpy:.2f} kJ/mol.")
    print("\n")
    return missing_enthalpy

# Perform analysis for both quantum dots
missing_H_large = analyze_quantum_dot("Largest", mass_fraction_oleate_large, dH_measured_large)
missing_H_small = analyze_quantum_dot("Smallest", mass_fraction_oleate_small, dH_measured_small)

# --- Conclusion ---
print("--- Conclusion ---")
print("The calculations show a very large discrepancy between the measured enthalpy and the value expected from simply combining InP dissolution and oleate protonation.")
print(f"For the largest QD, there is a missing endothermic contribution of {missing_H_large:.2f} kJ/mol.")
print(f"For the smallest QD, this missing endothermic contribution increases to {missing_H_small:.2f} kJ/mol.")
print("\nThis means there is a significant energy input required during dissolution that is not accounted for by the protonation of oleate alone (Choice A).")
print("This missing energy is larger for the smaller quantum dots, which have a higher surface-area-to-volume ratio and a greater proportion of ligands.")
print("Choice D provides the most logical explanation: disrupting the tightly packed shell of organic ligands requires a large amount of energy (endothermic), and this effect is more pronounced for smaller QDs which have a greater proportion of these ligands.")


# --- Final Output ---
# Restore stdout and print the captured output
sys.stdout = old_stdout
print(captured_output.getvalue())
print("<<<D>>>")