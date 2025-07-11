def analyze_heritability_statements():
    """
    Analyzes statements about heritability based on provided information.
    
    The script will:
    1. Define the parameters given in the problem (H^2 = 0.5).
    2. Model the components of variance (V_P, V_G, V_E, V_A, V_D, V_I).
    3. Evaluate each statement based on the principles of quantitative genetics.
    4. Print the step-by-step reasoning for each conclusion.
    """
    
    # --- 1. Setup based on the problem ---
    H2 = 0.5  # Broad-sense heritability (V_G / V_P)
    # For clarity, let's normalize total phenotypic variance to 1.0
    V_P = 1.0 
    V_G = H2 * V_P
    V_E = V_P - V_G

    print("--- Problem Analysis ---")
    print(f"Given: Broad-sense heritability (H^2) = V_G / V_P = {H2}")
    print(f"This means total Genetic Variance (V_G) is 50% of Phenotypic Variance (V_P).")
    print("Variance Equation: V_P = V_G + V_E")
    print(f"If we set V_P = {V_P}, then V_G = {V_G} and V_E = {V_E}.")
    print("\nGenetic Variance Equation: V_G = V_A (Additive) + V_D (Dominance) + V_I (Epistatic)")
    print("A standard Polygenic Score (PGS) estimates variance from additive effects, approaching h^2 = V_A / V_P.")
    print("-" * 28 + "\n")

    # --- 2. Evaluate Statement A ---
    print("--- Evaluating Statement A ---")
    print("Statement A: The PGS can not explain more than 50% of the variance.")
    print("The maximum variance a PGS can explain is h^2 = V_A / V_P.")
    print("The total genetic variance is H^2 = (V_A + V_D + V_I) / V_P = 0.5.")
    print("By definition, V_A must be less than or equal to V_G (since V_D and V_I cannot be negative).")
    print(f"Therefore, h^2 <= H^2, which means the explained variance must be <= {H2}.")
    print("Conclusion: Statement A is NECESSARILY TRUE.\n")

    # --- 3. Evaluate Statement C ---
    print("--- Evaluating Statement C ---")
    print("Statement C: The PGS will not approach 50% explained variance due to non-linear effects.")
    print("This statement implies that non-additive effects (V_D or V_I) are greater than zero.")
    # Example scenario with non-additive effects
    V_D_example = 0.1
    V_I_example = 0.05
    V_A_example = V_G - V_D_example - V_I_example
    h2_example = V_A_example / V_P
    print(f"If V_D or V_I > 0 (e.g., V_D={V_D_example}, V_I={V_I_example}), then V_A = {V_G} - {V_D_example} - {V_I_example} = {V_A_example:.2f}.")
    print(f"Then h^2 = {h2_example:.2f} / {V_P} = {h2_example:.2f}, which is less than H^2 ({H2}).")
    print("For a 'polygenic' trait, the existence of some non-additive effects is the standard biological assumption.")
    print("The question's phrasing focuses on the limitation of the 'linearly summing' method, which cannot capture these effects.")
    print("Conclusion: Statement C is considered TRUE in the context of a typical polygenic trait.\n")
    
    # --- 4. Evaluate Statement B & D ---
    print("--- Evaluating Statements B and D ---")
    print("Statement B (PGS will approach 50%) is false for the same reason C is true.")
    print("Statement D (Epigenetics limits h^2 < 0.5) is false because epigenetic effects are part of V_E and do not dictate the composition of V_G.\n")
    
    # --- 5. Final Conclusion ---
    print("--- Final Conclusion ---")
    print("Statement A is true by mathematical definition.")
    print("Statement C is true under the standard interpretation of a complex, polygenic trait.")
    print("Therefore, the option stating that A and C are correct is the best answer.")


if __name__ == "__main__":
    analyze_heritability_statements()