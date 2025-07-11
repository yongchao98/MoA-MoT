import math

def solve():
    """
    Analyzes the energetic contributions to the dissolution of InP quantum dots
    to determine the most logical explanation for the observed enthalpy changes.
    """
    # Print intro
    print("Analyzing the energetic contributions to the dissolution of InP quantum dots.")
    print("The goal is to determine if the increased amount of oleate ligand and its subsequent protonation can explain the observed increase in endothermicity for smaller quantum dots.")
    print("-" * 70)

    # Given constants and data
    mw_inp = 146  # g/mol
    mw_oleate = 281  # g/mol
    h_protonation_oleate = 7  # kJ/mol (endothermic)

    # Data for the largest quantum dot
    mass_frac_oleate_large = 0.46
    h_diss_large = 70  # kJ/mol of InP

    # Data for the smallest quantum dot
    mass_frac_oleate_small = 0.52
    h_diss_small = 120  # kJ/mol of InP

    # --- Calculations for the LARGEST Quantum Dot ---
    # Based on a 100g sample to find the molar ratio
    mass_inp_large = 100 * (1 - mass_frac_oleate_large)
    mass_oleate_large = 100 * mass_frac_oleate_large
    moles_inp_large = mass_inp_large / mw_inp
    moles_oleate_large = mass_oleate_large / mw_oleate
    # Moles of oleate per mole of InP
    ratio_large = moles_oleate_large / moles_inp_large
    # Enthalpy contribution from protonation per mole of InP
    h_prot_large = ratio_large * h_protonation_oleate

    # --- Calculations for the SMALLEST Quantum Dot ---
    # Based on a 100g sample to find the molar ratio
    mass_inp_small = 100 * (1 - mass_frac_oleate_small)
    mass_oleate_small = 100 * mass_frac_oleate_small
    moles_inp_small = mass_inp_small / mw_inp
    moles_oleate_small = mass_oleate_small / mw_oleate
    # Moles of oleate per mole of InP
    ratio_small = moles_oleate_small / moles_inp_small
    # Enthalpy contribution from protonation per mole of InP
    h_prot_small = ratio_small * h_protonation_oleate

    # --- Comparison of Changes ---
    observed_enthalpy_change = h_diss_small - h_diss_large
    protonation_enthalpy_change = h_prot_small - h_prot_large

    # --- Final Output ---
    print("--- Analysis Results ---")
    print(f"For the largest quantum dot:")
    print(f"The total enthalpy of dissolution is {h_diss_large} kJ/mol of InP.")
    print(f"The mole ratio of oleate to InP is {ratio_large:.2f}.")
    print(f"The calculated contribution from oleate protonation is {ratio_large:.2f} * {h_protonation_oleate} = {h_prot_large:.2f} kJ/mol of InP.")
    print("\n")

    print(f"For the smallest quantum dot:")
    print(f"The total enthalpy of dissolution is {h_diss_small} kJ/mol of InP.")
    print(f"The mole ratio of oleate to InP is {ratio_small:.2f}.")
    print(f"The calculated contribution from oleate protonation is {ratio_small:.2f} * {h_protonation_oleate} = {h_prot_small:.2f} kJ/mol of InP.")
    print("-" * 70)
    
    print("--- Comparison of Enthalpy Changes ---")
    print(f"The observed increase in dissolution enthalpy is: {h_diss_small} - {h_diss_large} = {observed_enthalpy_change:.2f} kJ/mol.")
    print(f"The calculated increase from oleate protonation is: {h_prot_small:.2f} - {h_prot_large:.2f} = {protonation_enthalpy_change:.2f} kJ/mol.")
    print("\n")
    print("Conclusion: The increased protonation energy only accounts for ~{protonation_enthalpy_change:.1f} kJ/mol of the observed {observed_enthalpy_change:.1f} kJ/mol change.")
    print("This demonstrates that while the protonation effect (Option A) is real, it is far too small to be the primary explanation.")
    print("Therefore, a much larger endothermic effect that scales with the surface area must be at play, such as the energy required to disrupt the packed ligand shell (Option D).")

solve()