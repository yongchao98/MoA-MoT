import sys
from io import StringIO

# Redirect stdout to a string buffer to capture the output
old_stdout = sys.stdout
sys.stdout = mystdout = StringIO()

def solve_enthalpy_puzzle():
    """
    Calculates the contribution of oleate protonation to the total enthalpy of dissolution
    for large and small quantum dots to test the hypothesis in answer choice A.
    """
    # Given values
    mw_inp = 146  # g/mol
    mw_oleate = 281  # g/mol
    h_prot_oleate = 7  # kJ/mol of oleate

    # Largest Quantum Dot Data
    mass_frac_oleate_large = 0.46
    mass_frac_inp_large = 1 - mass_frac_oleate_large
    h_diss_large = 70  # kJ/mol of InP

    # Smallest Quantum Dot Data
    mass_frac_oleate_small = 0.52
    mass_frac_inp_small = 1 - mass_frac_oleate_small
    h_diss_small = 120  # kJ/mol of InP

    # --- Calculations for the largest quantum dot ---
    # Consider a 100g sample to find the mole ratio
    moles_oleate_large = mass_frac_oleate_large * 100 / mw_oleate
    moles_inp_large = mass_frac_inp_large * 100 / mw_inp
    # Mole ratio of oleate to InP
    ratio_large = moles_oleate_large / moles_inp_large
    # Enthalpy contribution from oleate protonation (per mole of InP)
    h_contribution_large = ratio_large * h_prot_oleate

    # --- Calculations for the smallest quantum dot ---
    # Consider a 100g sample to find the mole ratio
    moles_oleate_small = mass_frac_oleate_small * 100 / mw_oleate
    moles_inp_small = mass_frac_inp_small * 100 / mw_inp
    # Mole ratio of oleate to InP
    ratio_small = moles_oleate_small / moles_inp_small
    # Enthalpy contribution from oleate protonation (per mole of InP)
    h_contribution_small = ratio_small * h_prot_oleate

    # --- Compare calculated change with observed change ---
    observed_h_change = h_diss_small - h_diss_large
    calculated_h_change_from_protonation = h_contribution_small - h_contribution_large

    # --- Print results ---
    print("Analysis of Enthalpy Contributions:")
    print("-" * 35)
    print(f"Observed total enthalpy change from large to small QD: {observed_h_change:.1f} kJ/mol")
    print("\n--- Contribution from Oleate Protonation ---")
    print(f"For the largest QD:")
    print(f"  Enthalpy from oleate protonation = {h_contribution_large:.2f} kJ/mol of InP")
    print(f"For the smallest QD:")
    print(f"  Enthalpy from oleate protonation = {h_contribution_small:.2f} kJ/mol of InP")
    print("\n" + "="*50)
    print(f"Change in enthalpy explained by oleate protonation alone: {h_contribution_small:.2f} - {h_contribution_large:.2f} = {calculated_h_change_from_protonation:.2f} kJ/mol")
    print("="*50 + "\n")
    print("Conclusion:")
    print(f"The protonation of oleate accounts for only {calculated_h_change_from_protonation:.2f} kJ/mol of the observed {observed_h_change:.1f} kJ/mol increase in endothermic enthalpy.")
    print("This is a very small fraction of the total change, so this effect alone cannot be the primary explanation.")

solve_enthalpy_puzzle()

# Restore stdout and print the captured output
sys.stdout = old_stdout
print(mystdout.getvalue())