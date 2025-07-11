def calculate_heritability(scenario, V_A, V_D, V_I, V_E):
    """
    Calculates and prints broad-sense (H^2) and narrow-sense (h^2) heritability.

    Args:
        scenario (str): A description of the scenario.
        V_A (float): Additive genetic variance.
        V_D (float): Dominance genetic variance.
        V_I (float): Epistatic (interaction) genetic variance.
        V_E (float): Environmental variance.
    """
    print(f"--- {scenario} ---")
    
    # Calculate total genetic and phenotypic variances
    V_G = V_A + V_D + V_I
    V_P = V_G + V_E
    
    # Avoid division by zero
    if V_P == 0:
        print("Total phenotypic variance is zero, heritability cannot be calculated.")
        return
        
    # Calculate broad-sense heritability (H^2)
    H_sq = V_G / V_P
    
    # Calculate narrow-sense heritability (h^2)
    h_sq = V_A / V_P
    
    print(f"Variance Components: V_A = {V_A}, V_D = {V_D}, V_I = {V_I}, V_E = {V_E}")
    print(f"Total Genetic Variance (V_G) = {V_A} + {V_D} + {V_I} = {V_G:.2f}")
    print(f"Total Phenotypic Variance (V_P) = V_G + V_E = {V_G:.2f} + {V_E} = {V_P:.2f}")
    print("\nCalculated Heritability:")
    print(f"Broad-sense Heritability (H^2) = V_G / V_P = {V_G:.2f} / {V_P:.2f} = {H_sq:.4f}")
    print(f"Narrow-sense Heritability (h^2) = V_A / V_P = {V_A} / {V_P:.2f} = {h_sq:.4f}\n")

    if abs(H_sq - h_sq) < 1e-9:
        print("Result: H^2 and h^2 are equal, as expected in a purely additive system.")
    else:
        print("Result: H^2 and h^2 are different due to non-additive genetic variance (V_D or V_I).")
    print("-" * 25)

# --- Scenario 1: Rabbit Experiment (Purely Additive Model) ---
# Given H^2 = 0.75 and the system is purely additive, so h^2 = 0.75.
# Let's assume V_P = 100. Then V_G = 75. Since it's all additive, V_A = 75.
# This leaves V_E = V_P - V_G = 25.
calculate_heritability(
    scenario="Rabbit Model (Purely Additive)",
    V_A=75.0, 
    V_D=0.0, 
    V_I=0.0, 
    V_E=25.0
)

# --- Scenario 2: Different Species (with Epistasis, as per Choice C) ---
# We introduce epistatic variance (V_I > 0) to see the effect.
# Let's keep V_A and V_E the same and add some V_I.
calculate_heritability(
    scenario="Other Species Model (with Epistasis)",
    V_A=75.0, 
    V_D=0.0, 
    V_I=15.0, # Introducing epistatic variance
    V_E=25.0
)