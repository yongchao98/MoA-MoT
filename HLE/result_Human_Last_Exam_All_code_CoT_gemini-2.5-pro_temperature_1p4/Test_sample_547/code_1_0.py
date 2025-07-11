import math

def analyze_dissolution_enthalpy():
    """
    Analyzes the enthalpy of dissolution for InP quantum dots of different sizes.
    """
    # Given constants
    mw_inp = 146  # g/mol
    mw_oleate = 281  # g/mol
    delta_h_diss_bulk_inp = -86  # kJ/mol
    delta_h_prot_oleate = 7  # kJ/mol

    # Properties of the largest quantum dot
    mass_frac_oleate_large = 0.46
    delta_h_diss_large_qd = 70  # kJ/mol of InP

    # Properties of the smallest quantum dot
    mass_frac_oleate_small = 0.52
    delta_h_diss_small_qd = 120  # kJ/mol of InP
    
    # --- Calculations for the largest quantum dot ---
    # Assume a 100g sample
    mass_oleate_large = 100 * mass_frac_oleate_large
    mass_inp_large = 100 - mass_oleate_large
    
    moles_oleate_large = mass_oleate_large / mw_oleate
    moles_inp_large = mass_inp_large / mw_inp
    
    # Molar ratio of oleate to InP
    ratio_large = moles_oleate_large / moles_inp_large
    
    # Enthalpy contribution from oleate protonation (per mole of InP)
    h_oleate_contribution_large = ratio_large * delta_h_prot_oleate
    
    # Calculate the excess surface enthalpy using the energy balance equation:
    # ΔH_total = ΔH_surface + ΔH(bulk_InP_diss) + n * ΔH(oleate_prot)
    # ΔH_surface = ΔH_total - ΔH(bulk_InP_diss) - n * ΔH(oleate_prot)
    h_surface_large = delta_h_diss_large_qd - delta_h_diss_bulk_inp - h_oleate_contribution_large

    # --- Calculations for the smallest quantum dot ---
    # Assume a 100g sample
    mass_oleate_small = 100 * mass_frac_oleate_small
    mass_inp_small = 100 - mass_oleate_small
    
    moles_oleate_small = mass_oleate_small / mw_oleate
    moles_inp_small = mass_inp_small / mw_inp
    
    # Molar ratio of oleate to InP
    ratio_small = moles_oleate_small / moles_inp_small
    
    # Enthalpy contribution from oleate protonation (per mole of InP)
    h_oleate_contribution_small = ratio_small * delta_h_prot_oleate
    
    # Calculate the excess surface enthalpy
    h_surface_small = delta_h_diss_small_qd - delta_h_diss_bulk_inp - h_oleate_contribution_small
    
    # Print the results
    print("--- Analysis of Dissolution Enthalpy ---\n")
    print("For the LARGEST quantum dot:")
    print(f"Molar ratio (oleate/InP): {ratio_large:.3f}")
    print(f"Enthalpy from oleate protonation: {h_oleate_contribution_large:.2f} kJ/mol of InP")
    print("\nEnergy Balance Equation (per mole of InP):")
    print("ΔH_total = ΔH_surface + ΔH_bulk_InP + ΔH_oleate")
    print(f"{delta_h_diss_large_qd} kJ/mol = ΔH_surface + ({delta_h_diss_bulk_inp} kJ/mol) + ({h_oleate_contribution_large:.2f} kJ/mol)")
    print(f"Calculated Excess Surface Enthalpy (ΔH_surface): {h_surface_large:.2f} kJ/mol of InP")
    
    print("\n" + "="*40 + "\n")
    
    print("For the SMALLEST quantum dot:")
    print(f"Molar ratio (oleate/InP): {ratio_small:.3f}")
    print(f"Enthalpy from oleate protonation: {h_oleate_contribution_small:.2f} kJ/mol of InP")
    print("\nEnergy Balance Equation (per mole of InP):")
    print("ΔH_total = ΔH_surface + ΔH_bulk_InP + ΔH_oleate")
    print(f"{delta_h_diss_small_qd} kJ/mol = ΔH_surface + ({delta_h_diss_bulk_inp} kJ/mol) + ({h_oleate_contribution_small:.2f} kJ/mol)")
    print(f"Calculated Excess Surface Enthalpy (ΔH_surface): {h_surface_small:.2f} kJ/mol of InP")
    
    print("\n" + "="*40 + "\n")
    
    # Comparison and Conclusion
    total_change = delta_h_diss_small_qd - delta_h_diss_large_qd
    change_due_to_oleate = h_oleate_contribution_small - h_oleate_contribution_large
    
    print("Conclusion:")
    print(f"The total change in dissolution enthalpy is {total_change} kJ/mol.")
    print(f"The change due to increased oleate protonation is only {change_due_to_oleate:.2f} kJ/mol.")
    print("This confirms that the protonation of oleate is not the primary reason for the large increase in endothermicity.")
    print("The primary reason is the large and size-dependent 'Excess Surface Enthalpy' term, which represents the energy needed to disrupt the stable, passivated surface of the quantum dot.")


analyze_dissolution_enthalpy()