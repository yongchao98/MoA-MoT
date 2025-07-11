def calculate_heritability():
    """
    Calculates and prints heritability measures for two scenarios
    to demonstrate the effect of non-additive genetic variance.
    """
    print("This script demonstrates why the presence of non-additive genetic variance (like epistasis) causes broad-sense (H^2) and narrow-sense (h^2) heritability to differ.")
    print("-" * 70)

    # --- Scenario 1: Rabbit Model (Purely Additive Genetics) ---
    print("\nScenario 1: A purely additive model (as described for the rabbits).")
    print("Here, Dominance (V_D) and Epistatic (V_I) variances are zero.")

    V_A_1 = 60  # Additive genetic variance
    V_D_1 = 0   # Dominance genetic variance
    V_I_1 = 0   # Epistatic genetic variance
    V_E_1 = 20  # Environmental variance

    # Calculate total genetic and phenotypic variances
    V_G_1 = V_A_1 + V_D_1 + V_I_1
    V_P_1 = V_G_1 + V_E_1

    # Calculate heritability
    H_sq_1 = V_G_1 / V_P_1
    h_sq_1 = V_A_1 / V_P_1

    print(f"\nEquations for Scenario 1:")
    print(f"Total Genetic Variance (V_G) = V_A + V_D + V_I = {V_A_1} + {V_D_1} + {V_I_1} = {V_G_1}")
    print(f"Total Phenotypic Variance (V_P) = V_G + V_E = {V_G_1} + {V_E_1} = {V_P_1}")
    print(f"Broad-Sense Heritability (H^2) = V_G / V_P = {V_G_1} / {V_P_1} = {H_sq_1:.2f}")
    print(f"Narrow-Sense Heritability (h^2) = V_A / V_P = {V_A_1} / {V_P_1} = {h_sq_1:.2f}")
    print("\nResult: In a purely additive model, H^2 is equal to h^2.")
    print("-" * 70)

    # --- Scenario 2: Model with Epistasis ---
    print("\nScenario 2: A model including epistatic interactions.")
    print("Here, we introduce a non-zero Epistatic variance (V_I).")

    V_A_2 = 60  # Additive genetic variance (kept the same)
    V_D_2 = 0   # Dominance genetic variance
    V_I_2 = 15  # Epistatic genetic variance (non-zero)
    V_E_2 = 20  # Environmental variance (kept the same)

    # Calculate total genetic and phenotypic variances
    V_G_2 = V_A_2 + V_D_2 + V_I_2
    V_P_2 = V_G_2 + V_E_2

    # Calculate heritability
    H_sq_2 = V_G_2 / V_P_2
    h_sq_2 = V_A_2 / V_P_2

    print(f"\nEquations for Scenario 2:")
    print(f"Total Genetic Variance (V_G) = V_A + V_D + V_I = {V_A_2} + {V_D_2} + {V_I_2} = {V_G_2}")
    print(f"Total Phenotypic Variance (V_P) = V_G + V_E = {V_G_2} + {V_E_2} = {V_P_2}")
    print(f"Broad-Sense Heritability (H^2) = V_G / V_P = {V_G_2} / {V_P_2} = {H_sq_2:.2f}")
    print(f"Narrow-Sense Heritability (h^2) = V_A / V_P = {V_A_2} / {V_P_2} = {h_sq_2:.2f}")
    print("\nResult: With epistasis, H^2 is now greater than h^2. This divergence is the key difference.")
    print("-" * 70)

    print("\nConclusion: The presence of epistatic interactions (a non-additive genetic component) causes a difference between H^2 and h^2. This supports answer C.")

if __name__ == "__main__":
    calculate_heritability()
<<<C>>>