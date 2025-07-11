import math

def predict_isomer_ratio():
    """
    Calculates the ratio of cis (A) and trans (B) products from the thermal
    electrocyclization of (2Z,4Z,6Z,8E)-deca-2,4,6,8-tetraene based on FMO theory.
    """

    # Step 1: Define constants for the calculation.
    # The energy difference (ΔΔG‡) between the transition states leading to the
    # trans and cis products. This is estimated from the steric strain of a
    # methyl group (A-value).
    delta_delta_g = 1.7  # units: kcal/mol

    # The ideal gas constant.
    R = 0.001987  # units: kcal/mol·K

    # The temperature for the thermal reaction, assuming standard room temperature.
    T = 298  # units: Kelvin

    # Step 2: Explain the theoretical basis for the calculation.
    print("--- Theoretical Analysis ---")
    print("1. The reaction is a thermal electrocyclization of an 8-pi-electron system.")
    print("2. FMO theory predicts a conrotatory ring closure.")
    print("3. The (8E)-isomer cyclizes to the trans-product (B) via a transition state with a Me-H steric clash.")
    print("4. The (8Z)-isomer (formed via thermal isomerization) cyclizes to the cis-product (A) via a transition state with a less severe H-H steric clash.")
    print("5. The product ratio is determined by the energy difference (ΔΔG‡) between these two transition states.")
    print("   Ratio A/B = exp(ΔΔG‡ / RT)")
    print("-" * 28)

    # Step 3: Perform the calculation.
    rt_value = R * T
    exponent = delta_delta_g / rt_value
    ratio_A_to_B = math.exp(exponent)

    # Step 4: Display the final equation with all the numerical values.
    print("--- Calculation ---")
    print(f"The final equation is: Ratio = exp(ΔΔG‡ / (R * T))")
    print(f"Plugging in the values:")
    print(f"Ratio = exp({delta_delta_g} / ({R} * {T}))")
    print(f"Ratio = exp({delta_delta_g} / {rt_value:.3f})")
    print(f"Ratio = exp({exponent:.3f})")
    print("-" * 28)

    # Step 5: Print the final result.
    print(f"The predicted ratio of A (cis) to B (trans) is {ratio_A_to_B:.1f} : 1.")
    
    # Step 6: Output the final answer in the specified format.
    print(f"\n<<<The ratio of A to B is {ratio_A_to_B:.1f}>>>")

predict_isomer_ratio()