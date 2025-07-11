def calculate_heritability():
    """
    Demonstrates the difference between broad-sense (H^2) and narrow-sense (h^2)
    heritability by modeling two scenarios.
    """

    print("--- Scenario 1: Rabbit Population (Purely Additive Genetics) ---")
    # Set variance components based on the problem. H^2 = 0.75.
    # Let's assume V_G = 75 and V_E = 25, so V_P = 100. H^2 = 75/100 = 0.75.
    # In this scenario, "entirely additive" means V_G = V_A.
    V_A_1 = 75
    V_D_1 = 0
    V_I_1 = 0
    V_E_1 = 25

    # Calculate total genetic and phenotypic variance
    V_G_1 = V_A_1 + V_D_1 + V_I_1
    V_P_1 = V_G_1 + V_E_1

    # Calculate heritability values
    H2_1 = V_G_1 / V_P_1
    h2_1 = V_A_1 / V_P_1

    print(f"Variance Components: V_A = {V_A_1}, V_D = {V_D_1}, V_I = {V_I_1}, V_E = {V_E_1}")
    print(f"Total Genetic Variance (V_G) = {V_A_1} + {V_D_1} + {V_I_1} = {V_G_1}")
    print(f"Total Phenotypic Variance (V_P) = {V_G_1} + {V_E_1} = {V_P_1}")
    print("\nCalculating heritability:")
    print(f"Broad-sense (H^2) = V_G / V_P = {V_G_1} / {V_P_1} = {H2_1:.2f}")
    print(f"Narrow-sense (h^2) = V_A / V_P = {V_A_1} / {V_P_1} = {h2_1:.2f}")
    print("Result: In a purely additive system, H^2 is equal to h^2.\n")


    print("--- Scenario 2: Other Species (with Epistatic Interactions, as per Choice C) ---")
    # Let's keep Additive and Environmental variance the same but add Epistatic variance.
    V_A_2 = 75
    V_D_2 = 0
    V_I_2 = 15  # Introduce epistatic variance
    V_E_2 = 25

    # Calculate total genetic and phenotypic variance
    V_G_2 = V_A_2 + V_D_2 + V_I_2
    V_P_2 = V_G_2 + V_E_2

    # Calculate heritability values
    H2_2 = V_G_2 / V_P_2
    h2_2 = V_A_2 / V_P_2

    print(f"Variance Components: V_A = {V_A_2}, V_D = {V_D_2}, V_I = {V_I_2}, V_E = {V_E_2}")
    print(f"Total Genetic Variance (V_G) = {V_A_2} + {V_D_2} + {V_I_2} = {V_G_2}")
    print(f"Total Phenotypic Variance (V_P) = {V_G_2} + {V_E_2} = {V_P_2}")
    print("\nCalculating heritability:")
    print(f"Broad-sense (H^2) = V_G / V_P = {V_G_2} / {V_P_2} = {H2_2:.3f}")
    print(f"Narrow-sense (h^2) = V_A / V_P = {V_A_2} / {V_P_2} = {h2_2:.3f}")
    print("Result: The presence of epistatic variance (V_I > 0) causes H^2 to be different from h^2.")

# Execute the function to see the output
calculate_heritability()