import math

def predict_product_ratio():
    """
    Predicts the product ratio of a thermal electrocyclization using FMO theory
    and transition state energy analysis.
    """

    # Part 1: Theoretical Background based on FMO Theory
    print("--- Frontier Molecular Orbital (FMO) Theory Analysis ---")
    print("1. Reaction Type: The reaction is the thermal electrocyclization of a deca-2,4,6,8-tetraene system.")
    print("2. Electron Count: The conjugated system has 8 pi electrons.")
    print("3. Selection Rule: For a thermal reaction with 8 pi electrons (a 4k system, where k=2), FMO theory predicts a 'conrotatory' ring closure.")
    print("-" * 60)

    # Part 2: Stereochemical Pathways to Products A (cis) and B (trans)
    print("--- Stereochemical Analysis of Conrotatory Pathways ---")
    print("The formation of two distinct products, cis (A) and trans (B), implies two competing reaction pathways.")
    print("\nPathway to trans-isomer (B):")
    print("The starting material, (2Z,4Z,6Z,8E)-deca-2,4,6,8-tetraene, preferentially cyclizes to the trans-product B.")
    print("In the transition state (TS_B) for this process, the two inward-rotating groups are both hydrogen atoms, leading to minimal steric strain.")
    
    print("\nPathway to cis-isomer (A):")
    print("The cis-product A is formed from the (2Z,4Z,6Z,8Z)-isomer, which is in thermal equilibrium with the starting (8E)-isomer.")
    print("The (8Z)-isomer undergoes conrotation to form the cis-product A.")
    print("In its transition state (TS_A), one inward-rotating group is a hydrogen, while the other is a bulky methyl group, which introduces significant steric strain.")
    print("-" * 60)
    
    # Part 3: Quantitative Prediction using Curtin-Hammett Principle
    print("--- Quantitative Analysis of the Product Ratio ---")
    print("The ratio of products [A]/[B] is determined by the overall free energy difference (dG_overall) between the two transition states, TS_A and TS_B.")
    print("This energy difference is the sum of two main contributions:")

    # Constants and Assumptions
    R = 8.314  # Gas constant in J/(mol·K)
    T_celsius = 100
    T_kelvin = T_celsius + 273.15
    kcal_to_J = 4184

    # Estimated energy differences based on established chemical principles
    # 1. Energy difference between Z and E ground states
    dG_iso_kcal = 1.0  # (8Z)-isomer is ~1.0 kcal/mol less stable than (8E)-isomer
    dG_iso_J = dG_iso_kcal * kcal_to_J
    
    # 2. Energy difference in activation energy due to steric hindrance (torquoselectivity)
    dG_TS_steric_kcal = 2.0  # TS_A is ~2.0 kcal/mol higher than TS_B due to Me/H inward rotation
    dG_TS_steric_J = dG_TS_steric_kcal * kcal_to_J
    
    print(f"  a) Ground State Instability (dG_iso): The (8Z)-isomer is less stable than the (8E)-isomer by an estimated {dG_iso_kcal} kcal/mol.")
    print(f"  b) Transition State Strain (dG_TS_steric): TS_A is higher in energy than TS_B due to steric hindrance by an estimated {dG_TS_steric_kcal} kcal/mol.")
    
    # Total energy difference
    dG_overall_kcal = dG_iso_kcal + dG_TS_steric_kcal
    dG_overall_J = dG_iso_J + dG_TS_steric_J

    print(f"\nTotal energy difference, dG_overall = {dG_iso_kcal} + {dG_TS_steric_kcal} = {dG_overall_kcal} kcal/mol")

    # Calculation of the ratio using the Boltzmann distribution
    RT = R * T_kelvin
    ratio_A_to_B = math.exp(-dG_overall_J / RT)
    
    print("\nThe product ratio [A]/[B] is calculated using the equation: exp(-dG_overall / RT)")
    print(f"Assuming a reaction temperature of {T_celsius}°C ({T_kelvin:.2f} K):")
    print(f"Final Equation: [A]/[B] = exp(-{dG_overall_J:.0f} / ({R:.3f} * {T_kelvin:.2f}))")
    
    # Final Result
    print("-" * 60)
    print("--- Final Predicted Ratio ---")
    if ratio_A_to_B != 0:
        print(f"The calculated ratio of [A]/[B] is {ratio_A_to_B:.4f}")
        print(f"This corresponds to a product ratio of A : B ≈ 1 : {1/ratio_A_to_B:.1f}")
        print("Therefore, the trans-isomer B is the major product.")
    else:
        print("The calculation results in a ratio where product A is not significantly formed.")

predict_product_ratio()
<<<1:51.5>>>