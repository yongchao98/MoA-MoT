import math

def predict_product_ratio():
    """
    Calculates the predicted product ratio of an electrocyclic reaction
    based on the difference in activation energies of competing pathways.
    """
    # Activation free energy for the pathway to cis-isomer A in kcal/mol
    # (from cyclization of the (...,8Z) isomer)
    dG_A_kcal = 25.7

    # Activation free energy for the pathway to trans-isomer B in kcal/mol
    # (from cyclization of the (...,8E) isomer)
    dG_B_kcal = 28.5

    # Reaction temperature in Celsius, from the experimental study
    T_celsius = 110.6

    # Gas constant in kcal/(mol*K)
    R_kcal = 0.001987

    # --- Calculations ---
    # Convert temperature to Kelvin
    T_kelvin = T_celsius + 273.15

    # Calculate the difference in activation energies
    delta_dG_kcal = dG_A_kcal - dG_B_kcal

    # Calculate RT in kcal/mol
    RT_kcal = R_kcal * T_kelvin

    # Calculate the ratio using the equation ratio = exp(-delta_dG / RT)
    ratio = math.exp(-delta_dG_kcal / RT_kcal)

    # --- Output the results ---
    print("Predicting the product ratio A/B for the thermal electrocyclization.")
    print("-" * 30)
    print("The final calculation is based on the equation: Ratio = exp(-(E_A - E_B) / (R * T))\n")

    print(f"Values used in the calculation:")
    print(f"Activation Energy for cis-product (A), E_A = {dG_A_kcal} kcal/mol")
    print(f"Activation Energy for trans-product (B), E_B = {dG_B_kcal} kcal/mol")
    print(f"Difference in activation energy (E_A - E_B) = {delta_dG_kcal:.2f} kcal/mol")
    print(f"Gas Constant, R = {R_kcal} kcal/(mol*K)")
    print(f"Temperature, T = {T_kelvin:.2f} K ({T_celsius}°C)")
    print("-" * 30)

    # Output the final equation with all the numbers
    print("Final equation with numbers:")
    print(f"Ratio = exp(-({dG_A_kcal} - {dG_B_kcal}) / ({R_kcal} * {T_kelvin:.2f}))")
    print(f"Ratio = exp(-({delta_dG_kcal:.2f}) / ({RT_kcal:.4f}))")
    print(f"Ratio = exp({-delta_dG_kcal:.2f} / {RT_kcal:.4f})")


    print(f"\nThe predicted ratio of cis-isomer (A) to trans-isomer (B) is: {ratio:.1f}")
    return ratio

# Run the calculation and print the final result
final_ratio = predict_product_ratio()

# The final answer in the required format
# <<<39.3>>>