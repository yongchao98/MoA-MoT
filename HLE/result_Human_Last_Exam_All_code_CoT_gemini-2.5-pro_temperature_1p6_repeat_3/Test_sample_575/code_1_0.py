def heritability_analysis():
    """
    Analyzes and explains the relationships between heritability types and
    the predictive power of a polygenic score (PGS).
    """
    H2 = 0.5  # Given broad-sense heritability

    # --- Introduction and Definitions ---
    print("Analyzing the limits of a Polygenic Score (PGS).")
    print(f"Given: Broad-sense heritability (H^2) = {H2}\n")

    print("Key Definitions:")
    print("1. Broad-Sense Heritability (H^2) = Total Genetic Variance (V_G) / Phenotypic Variance (V_P)")
    print("2. Narrow-Sense Heritability (h^2) = Additive Genetic Variance (V_A) / Phenotypic Variance (V_P)")
    print("3. Total Genetic Variance (V_G) = Additive (V_A) + Non-Additive Variance (V_D, V_I)")
    print("4. A standard PGS captures additive effects, so its maximum explained variance is h^2.\n")

    # --- The Logical Derivation ---
    print("--- Deriving the relationship between H^2 and h^2 ---")
    print("From the definition, we know that V_G = V_A + (Non-Additive Variance).")
    print("Since variance cannot be negative, the Non-Additive portion is >= 0.")
    print("Therefore, it must be true that: V_A <= V_G")
    print("Dividing by Phenotypic Variance (V_P): (V_A / V_P) <= (V_G / V_P)")
    print("This leads to the fundamental inequality: h^2 <= H^2\n")

    # --- Applying the problem's values ---
    print("--- Applying this to the problem ---")
    print(f"We are given H^2 = {H2}.")
    print(f"Based on our derivation, h^2 must be less than or equal to {H2}.")
    print(f"The maximum variance a PGS can explain is h^2.")
    print(f"Therefore, the maximum variance a PGS can explain is also less than or equal to {H2}.")
    print("\nConclusion: The statement 'The polygenic score can not explain more than 50% of the variance' is necessarily true.")

    # --- Illustrative Example to test other statements ---
    print("\n--- Testing other statements with a hypothetical example ---")
    V_P = 100
    V_G = H2 * V_P
    print(f"Let's assume total phenotypic variance V_P = {V_P}.")
    print(f"Then total genetic variance V_G = H^2 * V_P = {H2} * {V_P} = {V_G:.1f}\n")
    
    # Case 1: No non-additive effects
    print("Case 1: All genetic variance is additive (h^2 = H^2).")
    V_A_case1 = V_G
    h2_case1 = V_A_case1 / V_P
    print(f"  V_A = {V_A_case1:.1f}, Non-Additive Variance = 0.")
    print(f"  h^2 = V_A / V_P = {V_A_case1:.1f} / {V_P:.1f} = {h2_case1}")
    print("  In this case, a perfect PGS *would* approach explaining 50% of the variance. This shows Statement C is not necessarily true.\n")

    # Case 2: With non-additive effects
    print("Case 2: Some genetic variance is non-additive (h^2 < H^2).")
    V_A_case2 = 30.0
    V_non_additive = V_G - V_A_case2
    h2_case2 = V_A_case2 / V_P
    print(f"  Let V_A = {V_A_case2:.1f}, then Non-Additive Variance = {V_non_additive:.1f}.")
    print(f"  h^2 = V_A / V_P = {V_A_case2:.1f} / {V_P:.1f} = {h2_case2}")
    print("  In this case, a perfect PGS would approach explaining 30% of the variance, not 50%. This shows Statement B is not necessarily true.")


heritability_analysis()
<<<A>>>