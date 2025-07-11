import math

def predict_product_ratio():
    """
    This script predicts the product ratio of a thermal 8-pi electron electrocyclization
    based on the energy difference between the competing transition states.
    """
    
    # Plan:
    # 1. State the principle: The ratio of products A (cis) and B (trans) is
    #    determined by the difference in activation energies (dG_ddagger) of the
    #    two pathways, following the Boltzmann distribution: Ratio(A/B) = exp(-dG_ddagger / RT).
    # 2. Estimate dG_ddagger: The transition state (TS) leading to cis-isomer A
    #    has a severe Me-Me steric interaction, while the TS to trans-isomer B
    #    has a smaller Me-H interaction. The energy difference is estimated
    #    from standard allylic strain values.
    #    dG_ddagger = G_TS(cis) - G_TS(trans) approx 3.6 kcal/mol.
    # 3. Calculate the ratio. Product B (trans) is the major product as its
    #    transition state is lower in energy. We will calculate the B/A ratio.

    # --- Parameters ---
    # Energy difference between the transition states in kcal/mol
    delta_g_kcal_mol = 3.6
    
    # Gas constant in J/(mol*K)
    R = 8.314
    
    # Temperature in Kelvin (standard room temperature)
    T = 298
    
    # Conversion factor from kcal to J
    KCAL_TO_J = 4184

    # --- Calculation ---
    # Convert energy difference to J/mol
    delta_g_j_mol = delta_g_kcal_mol * KCAL_TO_J
    
    # The activation energy for the cis product is higher, so it's the minor product.
    # We calculate the ratio of major (B) to minor (A).
    # Ratio B/A = exp(- (G_TS(B) - G_TS(A)) / RT) = exp( (G_TS(A) - G_TS(B)) / RT)
    # Ratio B/A = exp(delta_g_j_mol / (R * T))
    ratio_B_over_A = math.exp(delta_g_j_mol / (R * T))

    # --- Output ---
    print("Prediction of the product ratio for the electrocyclization of (2Z,4Z,6Z,8E)-deca-2,4,6,8-tetraene.")
    print("-" * 80)
    print("The ratio of trans-isomer B to cis-isomer A is determined by the equation:")
    print(f"Ratio(B/A) = exp(ΔG‡ / (R * T))")
    print("\nUsing the following values:")
    print(f"ΔG‡ (energy difference between TS_cis and TS_trans) = {delta_g_kcal_mol} kcal/mol = {delta_g_j_mol} J/mol")
    print(f"R (gas constant) = {R} J/mol·K")
    print(f"T (temperature) = {T} K")

    print("\nThe calculation is:")
    print(f"Ratio(B/A) = exp({delta_g_j_mol:.1f} / ({R} * {T}))")
    
    print("\nResult:")
    print(f"The predicted ratio of trans-isomer B to cis-isomer A is approximately {ratio_B_over_A:.0f} : 1.")
    print(f"This means the ratio of A : B is 1 : {ratio_B_over_A:.0f}.")
    print("-" * 80)
    
    # Return the numerical answer for the ratio B/A
    return ratio_B_over_A

# Execute the function and print the final result.
final_ratio = predict_product_ratio()
<<<437>>>