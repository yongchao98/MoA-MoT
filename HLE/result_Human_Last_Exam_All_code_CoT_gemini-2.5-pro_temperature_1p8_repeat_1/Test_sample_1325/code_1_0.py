def calculate_heritability():
    """
    Demonstrates the difference between broad-sense (H^2) and narrow-sense (h^2)
    heritability by modeling two scenarios based on the problem.
    """
    
    # Case 1: Rabbit Experiment (as described in the problem)
    # Given H^2 = 0.75 and genetic variance is "entirely additive".
    # This implies Dominance (V_D) and Epistatic (V_I) variances are 0.
    # We can assume a Phenotypic Variance (V_P) of 100 for a simple example.
    print("--- Case 1: Rabbit Experiment (Entirely Additive) ---")
    H2_case1 = 0.75
    V_P_case1 = 100.0
    V_A_case1 = H2_case1 * V_P_case1
    V_D_case1 = 0.0
    V_I_case1 = 0.0
    V_G_case1 = V_A_case1 + V_D_case1 + V_I_case1
    V_E_case1 = V_P_case1 - V_G_case1

    h2_case1 = V_A_case1 / V_P_case1
    
    print(f"Given Components: V_A={V_A_case1}, V_D={V_D_case1}, V_I={V_I_case1}, V_E={V_E_case1}")
    print(f"Total Phenotypic Variance V_P = {V_P_case1}")
    print("\nBroad-Sense Heritability (H^2) = (V_A + V_D + V_I) / V_P")
    print(f"H^2 = ({V_A_case1} + {V_D_case1} + {V_I_case1}) / {V_P_case1} = {V_G_case1 / V_P_case1:.2f}")
    
    print("\nNarrow-Sense Heritability (h^2) = V_A / V_P")
    print(f"h^2 = {V_A_case1} / {V_P_case1} = {h2_case1:.2f}")
    print("\nResult: When genetic variance is purely additive, H^2 equals h^2.")
    print("-" * 50)
    
    # Case 2: Other Species (with Dominance Variance introduced)
    # This scenario demonstrates the principle in answer choice E.
    # We add dominance variance (V_D > 0).
    print("--- Case 2: Other Species (With Dominance Variance) ---")
    V_A_case2 = V_A_case1  # Keep additive variance the same for comparison
    V_D_case2 = 15.0       # Introduce dominance variance
    V_I_case2 = 0.0        # Assume epistasis is negligible ("additive model")
    V_E_case2 = V_E_case1  # Keep environmental variance the same

    V_G_case2 = V_A_case2 + V_D_case2 + V_I_case2
    V_P_case2 = V_G_case2 + V_E_case2
    
    H2_case2 = V_G_case2 / V_P_case2
    h2_case2 = V_A_case2 / V_P_case2

    print(f"Given Components: V_A={V_A_case2}, V_D={V_D_case2}, V_I={V_I_case2}, V_E={V_E_case2}")
    print(f"New Total Phenotypic Variance V_P = {V_P_case2}")
    
    print("\nBroad-Sense Heritability (H^2) = (V_A + V_D + V_I) / V_P")
    print(f"H^2 = ({V_A_case2} + {V_D_case2} + {V_I_case2}) / {V_P_case2} = {H2_case2:.4f}")

    print("\nNarrow-Sense Heritability (h^2) = V_A / V_P")
    print(f"h^2 = {V_A_case2} / {V_P_case2} = {h2_case2:.4f}")
    print("\nResult: The presence of dominance variance causes H^2 to be greater than h^2.")

# Execute the function to print the calculations
calculate_heritability()