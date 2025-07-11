import sys

def solve():
    """
    Calculates the enthalpy contribution from oleate protonation to test if it
    explains the observed trend in quantum dot dissolution enthalpy.
    """
    # Given constants
    MW_InP = 146  # g/mol
    MW_oleate = 281  # g/mol
    H_protonation_oleate = 7  # kJ/mol
    H_diss_large_QD = 70  # kJ/mol of InP
    H_diss_small_QD = 120  # kJ/mol of InP
    mass_frac_oleate_large = 0.46
    mass_frac_oleate_small = 0.52

    # --- Analysis for the largest quantum dot ---
    mass_frac_InP_large = 1 - mass_frac_oleate_large
    # To find the molar ratio, consider a 100g sample.
    # The ratio of masses is (100 * mass_frac_oleate_large) to (100 * mass_frac_InP_large).
    # The molar ratio is (mass_oleate / MW_oleate) / (mass_InP / MW_InP).
    molar_ratio_large = (mass_frac_oleate_large / MW_oleate) / (mass_frac_InP_large / MW_InP)
    h_contribution_large = molar_ratio_large * H_protonation_oleate

    print("--- Analysis of Oleate Protonation Enthalpy ---")
    print("\nFor the largest quantum dot (46% oleate):")
    sys.stdout.write(f"The molar ratio of oleate to InP is ({mass_frac_oleate_large:.2f} / {MW_oleate}) / ({mass_frac_InP_large:.2f} / {MW_InP})")
    print(f" = {molar_ratio_large:.3f} mol oleate per mol InP.")
    sys.stdout.write(f"The resulting enthalpy contribution from protonation is {molar_ratio_large:.3f} * {H_protonation_oleate}")
    print(f" = {h_contribution_large:.2f} kJ per mol of InP.")

    # --- Analysis for the smallest quantum dot ---
    mass_frac_InP_small = 1 - mass_frac_oleate_small
    molar_ratio_small = (mass_frac_oleate_small / MW_oleate) / (mass_frac_InP_small / MW_InP)
    h_contribution_small = molar_ratio_small * H_protonation_oleate
    
    print("\nFor the smallest quantum dot (52% oleate):")
    sys.stdout.write(f"The molar ratio of oleate to InP is ({mass_frac_oleate_small:.2f} / {MW_oleate}) / ({mass_frac_InP_small:.2f} / {MW_InP})")
    print(f" = {molar_ratio_small:.3f} mol oleate per mol InP.")
    sys.stdout.write(f"The resulting enthalpy contribution from protonation is {molar_ratio_small:.3f} * {H_protonation_oleate}")
    print(f" = {h_contribution_small:.2f} kJ per mol of InP.")

    # --- Comparison and Conclusion ---
    change_in_h_contribution = h_contribution_small - h_contribution_large
    observed_change_in_h = H_diss_small_QD - H_diss_large_QD
    
    print("\n--- Comparison ---")
    sys.stdout.write(f"The change in enthalpy due to oleate protonation is {h_contribution_small:.2f} - {h_contribution_large:.2f}")
    print(f" = {change_in_h_contribution:.2f} kJ/mol.")
    sys.stdout.write(f"The total observed change in dissolution enthalpy is {H_diss_small_QD} - {H_diss_large_QD}")
    print(f" = {observed_change_in_h:.2f} kJ/mol.")
    
    print("\nConclusion:")
    print("The enthalpy change from oleate protonation accounts for only a very small fraction")
    print("of the total observed change, so it cannot be the main explanation.")


solve()