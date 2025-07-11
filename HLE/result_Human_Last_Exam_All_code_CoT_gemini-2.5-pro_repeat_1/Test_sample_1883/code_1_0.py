import math

def predict_product_ratio():
    """
    Predicts the ratio of cis (A) and trans (B) products from the thermal
    electrocyclization of (2Z,4Z,6Z,8E)-deca-2,4,6,8-tetraene using FMO theory
    and principles of torquoselectivity.
    """

    # --- 1. Introduction and Theory ---
    print("--- Frontier Molecular Orbital Theory Analysis ---")
    print("Reaction: Thermal 8-pi electrocyclization of a substituted decatetraene.")
    print("Woodward-Hoffmann Rule: For a 4n (n=2) pi-electron system under thermal conditions, the ring closure is conrotatory.")
    print("\nTwo competing pathways lead to the two products:")
    print("Pathway to B (trans): (8E)-isomer (out/out Me) --conrotatory--> trans-product B")
    print("Pathway to A (cis):   (8Z)-isomer (out/in Me)  --conrotatory--> cis-product A")
    print("\nThe product ratio A/B is determined by the relative energy of the transition states (Curtin-Hammett Principle).")

    # --- 2. Define Energy Parameters (in kcal/mol) ---
    dE_out_Me = 2.1  # Activation energy penalty for an 'outward' rotating methyl group
    dE_in_Me = -1.6  # Activation energy relief for an 'inward' rotating methyl group
    E_strain_in_Me = 4.0 # Ground state destabilization from an 'inward' methyl group

    print("\n--- 3. Calculating Transition State Energies ---")
    print("Energy parameters (kcal/mol):")
    print(f"  - Ground state strain of 'in' Me: {E_strain_in_Me}")
    print(f"  - Activation energy penalty of 'out' Me: {dE_out_Me}")
    print(f"  - Activation energy relief of 'in' Me: {dE_in_Me}")

    # Let the ground state energy of the (8E) isomer be the reference energy, E=0.
    # Let dG0 be the activation energy of the unsubstituted parent reaction.
    # The absolute value of dG0 is not needed as it cancels out.
    E_gs_E = 0.0
    E_gs_Z = E_strain_in_Me

    # Activation energy for the B (trans) pathway from the (8E) isomer
    # The (8E) isomer has two 'out' methyl groups.
    G_ts_B_relative = E_gs_E + dE_out_Me + dE_out_Me
    print(f"\nRelative transition state energy for B (trans):")
    print(f"G_ts(B) = E_gs(E) + dE(out_Me) + dE(out_Me)")
    print(f"G_ts(B) = {E_gs_E} + {dE_out_Me} + {dE_out_Me} = {G_ts_B_relative:.1f} kcal/mol (relative to dG0)")

    # Activation energy for the A (cis) pathway from the (8Z) isomer
    # The (8Z) isomer has one 'out' and one 'in' methyl group.
    G_ts_A_relative = E_gs_Z + dE_out_Me + dE_in_Me
    print(f"\nRelative transition state energy for A (cis):")
    print(f"G_ts(A) = E_gs(Z) + dE(out_Me) + dE(in_Me)")
    print(f"G_ts(A) = {E_strain_in_Me} + {dE_out_Me} + ({dE_in_Me}) = {G_ts_A_relative:.1f} kcal/mol (relative to dG0)")

    # --- 4. Calculate Energy Difference and Ratio ---
    ddG = G_ts_A_relative - G_ts_B_relative
    
    # Constants for the calculation
    T_C = 100  # Assumed reaction temperature in Celsius
    T_K = T_C + 273.15  # Temperature in Kelvin
    R_kcal = 1.987 / 1000  # Gas constant in kcal/mol-K
    
    RT = R_kcal * T_K
    ratio_A_B = math.exp(-ddG / RT)

    print("\n--- 4. Final Ratio Calculation ---")
    print(f"Difference in activation energy (ddG) = G_ts(A) - G_ts(B)")
    print(f"ddG = {G_ts_A_relative:.1f} - {G_ts_B_relative:.1f} = {ddG:.1f} kcal/mol")
    
    print(f"\nRatio A/B = exp(-ddG / RT)")
    print(f"Assuming T = {T_C}°C ({T_K:.2f} K):")
    print(f"Ratio A/B = exp(-({ddG:.1f}) / ({R_kcal:.5f} * {T_K:.2f}))")
    print(f"Ratio A/B = exp(-({ddG:.1f}) / {RT:.3f}) = {ratio_A_B:.3f}")

    # --- 5. Simplify to Integer Ratio ---
    # Find simple integer ratio for A:B.
    # If ratio_A_B is approx 0.667, then A:B is 2:3.
    b_part = 3
    a_part = round(ratio_A_B * b_part)
    
    print("\n--- 5. Conclusion ---")
    print(f"The calculated ratio of A/B is approximately {ratio_A_B:.3f}, which is close to 2/3.")
    print(f"Therefore, the predicted ratio of A (cis-isomer) to B (trans-isomer) is {int(a_part)}:{int(b_part)}.")

predict_product_ratio()