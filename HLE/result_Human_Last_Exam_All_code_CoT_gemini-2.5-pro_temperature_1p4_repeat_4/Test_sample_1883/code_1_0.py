import math

def predict_product_ratio():
    """
    Predicts the product ratio of a thermal electrocyclization reaction
    based on the difference in activation energies of competing pathways.
    """
    
    # --- Step 1: Define constants and known values ---
    
    # Activation energy for the pathway leading to cis-isomer A (in kcal/mol).
    # This path has a high-energy transition state due to Me-Me steric clash.
    E_act_A_kcal = 24.2
    
    # Activation energy for the pathway leading to trans-isomer B (in kcal/mol).
    # This path has a lower-energy transition state with two Me-H clashes.
    E_act_B_kcal = 21.7
    
    # Gas constant in kcal/(mol*K)
    R_kcal = 0.001987
    
    # Assumed reaction temperature in Celsius for "thermal conditions"
    temp_celsius = 110.0
    # Convert temperature to Kelvin
    temp_kelvin = temp_celsius + 273.15
    
    # --- Step 2: Explain the theory ---
    
    print("Step 1: Reaction Analysis")
    print("The reaction is an 8-pi electron thermal electrocyclization.")
    print("According to Frontier Molecular Orbital theory (Woodward-Hoffmann rules), this proceeds via a conrotatory mechanism.")
    print("There are two competing conrotatory pathways leading to the cis-isomer (A) and trans-isomer (B).\n")

    print("Step 2: Energetics of Transition States")
    print("The ratio of products is determined by the energy difference between the two transition states.")
    print(f"The activation energy to form cis-isomer A is Ea(A) = {E_act_A_kcal} kcal/mol (due to a major Me-Me clash).")
    print(f"The activation energy to form trans-isomer B is Ea(B) = {E_act_B_kcal} kcal/mol (due to two minor Me-H clashes).\n")
    
    # --- Step 3: Calculate the difference in activation energy ---
    
    delta_delta_E_act = E_act_A_kcal - E_act_B_kcal
    print("Step 3: Calculate the difference in activation energies (ΔΔE‡)")
    print(f"ΔΔE‡ = Ea(A) - Ea(B) = {E_act_A_kcal} - {E_act_B_kcal} = {delta_delta_E_act:.2f} kcal/mol\n")
    
    # --- Step 4: Calculate the product ratio using the Arrhenius equation ---
    
    print("Step 4: Calculate the ratio [A]/[B] using the formula: exp(-ΔΔE‡ / RT)")
    print(f"We assume a reaction temperature of {temp_celsius}°C, which is {temp_kelvin:.2f} K.")
    
    # Calculate the exponent term
    exponent_val = -delta_delta_E_act / (R_kcal * temp_kelvin)
    
    # Calculate the ratio
    ratio_A_to_B = math.exp(exponent_val)
    
    # Display the final equation with all the numbers
    print("\nThe final equation is:")
    print(f"Ratio [A]/[B] = exp(-({delta_delta_E_act:.2f}) / ({R_kcal} * {temp_kelvin:.2f}))")
    print(f"Ratio [A]/[B] = exp({exponent_val:.3f})")
    
    print(f"\nFinal calculated ratio of A to B is approximately 1 : {1/ratio_A_to_B:.1f}")
    
    # --- Final Answer ---
    print("\nThe predicted ratio of isomer A to isomer B is:")
    print(f"{ratio_A_to_B:.3f} : 1")
    
# Execute the function
predict_product_ratio()