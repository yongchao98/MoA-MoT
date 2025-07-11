import math

def predict_product_ratio():
    """
    Predicts the ratio of cis (A) and trans (B) products from the thermal
    electrocyclization of (2Z,4Z,6Z,8E)-deca-2,4,6,8-tetraene.
    """
    # Step 1: Explain the governing principle from FMO theory.
    print("Step 1: Determine the stereochemical mode of reaction.")
    print("The reaction is a thermal electrocyclization of an 8-pi electron system.")
    print("According to Woodward-Hoffmann rules, this is a 4n system (n=2), so the ring closure is conrotatory.\n")

    # Step 2: Explain the origin of the two products.
    print("Step 2: Identify the reaction pathways.")
    print("Two products (cis-A and trans-B) are formed because the reaction proceeds from two different ground-state conformers in equilibrium:")
    print("  - A low-energy 'helical' conformer cyclizes to give the trans-product B.")
    print("  - A higher-energy 'kinked' conformer cyclizes to give the cis-product A.\n")

    # Step 3: Relate product ratio to energy difference.
    print("Step 3: Determine the energy difference between pathways.")
    print("The transition states for both pathways involve a similar H-H steric clash.")
    print("Therefore, the difference in activation energy (ΔE) is approximately the energy difference between the two ground-state conformers.\n")

    # Step 4: Define constants and the estimated energy difference.
    print("Step 4: Set up the calculation.")
    # The energy difference is based on a standard value for steric strain.
    delta_E_kcal_mol = 0.9  # Assumed energy difference in kcal/mol
    # Conversion factor from kcal to J
    KCAL_TO_J = 4184
    delta_E_J_mol = delta_E_kcal_mol * KCAL_TO_J

    R = 8.314  # Gas constant in J/(mol*K)
    T = 298    # Standard temperature in K (25 °C)
    
    print(f"The energy difference ΔE is estimated as {delta_E_kcal_mol} kcal/mol ({delta_E_J_mol:.0f} J/mol).")
    print(f"Using R = {R} J/(mol*K) and T = {T} K.\n")

    # Step 5: Calculate the product ratio using the Boltzmann distribution.
    print("Step 5: Calculate the ratio [A]/[B].")
    print("The ratio is given by the Boltzmann distribution: Ratio = exp(-ΔE / RT)")

    exponent = -delta_E_J_mol / (R * T)
    ratio_A_to_B = math.exp(exponent)

    # Print the equation with the numbers plugged in.
    print(f"Ratio [A]/[B] = exp(-{delta_E_J_mol:.0f} / ({R:.3f} * {T}))")
    print(f"Ratio [A]/[B] = exp({exponent:.3f})")
    print(f"Ratio [A]/[B] = {ratio_A_to_B:.3f}\n")
    
    print(f"This means the trans-isomer (B) is the major product.")
    print(f"The predicted ratio of cis-isomer A to trans-isomer B is {ratio_A_to_B:.3f} : 1, or approximately 1 : {1/ratio_A_to_B:.1f}.")

predict_product_ratio()