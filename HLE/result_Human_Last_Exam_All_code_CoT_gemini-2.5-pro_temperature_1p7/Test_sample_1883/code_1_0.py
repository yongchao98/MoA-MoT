import math

def predict_product_ratio():
    """
    This function predicts the product ratio of a thermal electrocyclic reaction
    based on the principles of torquoselectivity and the Arrhenius equation.
    """

    # --- Step 1: Define Constants and Assumptions ---

    # Activation energy difference (delta_Ea) between the two pathways.
    # This value represents the steric penalty for rotating a methyl group inward.
    # We assume a literature-based value of 2.3 kcal/mol.
    delta_Ea_kcal_per_mol = 2.3

    # Assumed temperature (T) for the thermal reaction in Celsius.
    # 100 C is a common temperature for such reactions.
    T_celsius = 100.0

    # Ideal Gas Constant (R) in J/(mol*K).
    R_joules_per_mol_K = 8.314

    # Conversion factor from kcal to Joules.
    joules_per_kcal = 4184

    # --- Step 2: Convert Units ---

    # Convert delta_Ea from kcal/mol to J/mol.
    delta_Ea_joules_per_mol = delta_Ea_kcal_per_mol * joules_per_kcal

    # Convert temperature from Celsius to Kelvin.
    T_kelvin = T_celsius + 273.15

    # --- Step 3: Calculate the Ratio using the Arrhenius relationship ---
    # The ratio of products A/B = exp(delta_Ea / (R * T))
    # where delta_Ea = E_a(B) - E_a(A) is positive since A is the major product.

    exponent = delta_Ea_joules_per_mol / (R_joules_per_mol_K * T_kelvin)
    ratio_A_to_B = math.exp(exponent)

    # --- Step 4: Print the Explanation and Result ---

    print("Prediction of the product ratio for the thermal electrocyclization of (2Z,4Z,6Z,8E)-deca-2,4,6,8-tetraene.")
    print("The cis-isomer (A) is predicted to be the major product due to lower steric strain in its transition state.")
    print("\nThe product ratio (A/B) is calculated using the equation: Ratio = exp(ΔEa / (R * T))\n")

    print("The values used in the equation are:")
    print(f"Activation Energy Difference (ΔEa) = {delta_Ea_kcal_per_mol} kcal/mol or {delta_Ea_joules_per_mol:.1f} J/mol")
    print(f"Gas Constant (R) = {R_joules_per_mol_K} J/mol·K")
    print(f"Assumed Temperature (T) = {T_celsius}°C or {T_kelvin} K\n")

    print(f"The calculated ratio of cis-isomer (A) to trans-isomer (B) is approximately:")
    print(f"{ratio_A_to_B:.1f} : 1")


# Execute the function to get the prediction.
predict_product_ratio()
<<<22.3>>>