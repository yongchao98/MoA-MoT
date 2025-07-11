import math

def predict_product_ratio():
    """
    This function predicts the product ratio for the thermal electrocyclization of
    (2Z,4Z,6Z,8E)-deca-2,4,6,8-tetraene using Frontier Molecular Orbital theory
    and the principle of torquoselectivity.
    """

    # Step 1 & 2: Define the reaction and governing rule
    print("Step 1: Reaction Analysis")
    print("The reaction is a thermal electrocyclization of an 8-pi electron system.")
    print("According to the Woodward-Hoffmann rules, a thermal 8-pi (4n) electron reaction proceeds via a conrotatory mechanism.\n")

    # Step 3: Explain the formation of two products
    print("Step 2: Origin of the Two Products (A and B)")
    print("The starting material (2Z,4Z,6Z,8E)-tetraene cyclizes via conrotation to form the trans-isomer (Product B).")
    print("To form the cis-isomer (Product A), the starting material must first isomerize to the (2Z,4Z,6Z,8Z)-tetraene, which then undergoes conrotatory cyclization.")
    print("Product A (cis) is formed from the (2Z,8Z) isomer.")
    print("Product B (trans) is formed from the (2Z,8E) isomer.\n")

    # Step 4: Use torquoselectivity to compare transition state energies
    print("Step 3: FMO Theory and Torquoselectivity")
    print("The ratio of products A:B is determined by the relative energy of the transition states.")
    print("Torquoselectivity predicts that the pathway with the least steric hindrance is favored.")
    print("  - Pathway to A (cis): The most favorable transition state involves rotating two small Hydrogen atoms inward.")
    print("  - Pathway to B (trans): The transition state requires rotating one bulky Methyl group and one Hydrogen atom inward.")
    print("The transition state for product A is lower in energy because it is sterically less hindered.\n")

    # Step 5: Define constants and calculate the ratio
    # Delta_Ea: Energy difference between rotating a Me group vs. an H group inward (kcal/mol)
    # This is the difference in activation energies between the two pathways.
    delta_Ea = 1.4  # kcal/mol
    # R: Gas constant in kcal/mol-K
    R = 1.987 / 1000
    # T: Temperature in Kelvin (assuming room temperature)
    T = 298.0

    print("Step 4: Calculating the Ratio")
    print(f"The difference in activation energy (ΔEa) between the two pathways is estimated to be {delta_Ea} kcal/mol.")
    print(f"The reaction temperature (T) is assumed to be {T} K.")
    print(f"The gas constant (R) is {R:.6f} kcal/mol·K.")

    # Calculate the ratio using the Arrhenius relationship: Ratio = exp(ΔEa / RT)
    # We use Delta_Ea directly as it represents Ea(disfavored) - Ea(favored)
    ratio = math.exp(delta_Ea / (R * T))

    print("\nThe ratio of A:B is calculated as: exp(ΔEa / (R * T))")
    print(f"Ratio(A/B) = exp({delta_Ea} / ({R:.6f} * {T}))")
    print(f"Ratio(A/B) = exp({delta_Ea} / {R*T:.4f})")
    print(f"Ratio(A/B) = {ratio:.2f}\n")
    print(f"The predicted ratio of cis-isomer A to trans-isomer B is approximately {ratio:.1f} : 1.")
    
    # Final answer in the specified format
    print(f"\nFinal Answer in numeric format:")
    print(f'<<<{ratio:.1f}>>>')

# Run the prediction
predict_product_ratio()