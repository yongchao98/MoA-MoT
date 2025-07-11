def analyze_heritability_question():
    """
    Analyzes a quantitative genetics question to determine the correct statement.
    """
    # Step 1: Define the given information and key concepts.
    # The problem states the broad-sense heritability (H^2) is 0.5.
    # H^2 = Vg / Vp (Total Genetic Variance / Total Phenotypic Variance)
    H_squared = 0.5

    print("--- Problem Analysis ---")
    print(f"Given information: Broad-sense heritability (H^2) = {H_squared}")
    print("\nKey Concepts in Quantitative Genetics:")
    print("1. Broad-sense heritability (H^2): The proportion of phenotypic variance explained by ALL genetic factors (additive, dominance, epistasis).")
    print("2. Narrow-sense heritability (h^2): The proportion of phenotypic variance explained by ONLY ADDITIVE genetic factors.")
    print("3. Polygenic Score (PGS): A score built from GWAS data. Standard GWAS models and PGS are linear, so they primarily capture additive genetic effects.")

    # Step 2: Establish the mathematical relationships.
    # The total genetic variance (Vg) is the sum of additive (Va) and non-additive variance (V_non-additive).
    # Vg = Va + V_non-additive
    # Since V_non-additive >= 0, it is always true that Va <= Vg.
    # Dividing by total phenotypic variance (Vp), we get: Va/Vp <= Vg/Vp, which means h^2 <= H^2.
    # The maximum variance a standard PGS can theoretically explain (R^2_PGS) is the narrow-sense heritability (h^2).
    # So, R^2_PGS <= h^2.

    # Step 3: Combine relationships to find the upper bound for PGS performance.
    # From the above, we can chain the inequalities together.
    print("\nDeriving the Master Inequality:")
    print("Variance explained by PGS (R^2_PGS) is limited by narrow-sense heritability (h^2).")
    print("Narrow-sense heritability (h^2) is limited by broad-sense heritability (H^2).")
    print("\nThis gives us the final equation:")
    # Per instructions, printing each number in the final equation.
    print(f"R^2_PGS <= h^2 <= H^2")
    print(f"Substituting the given value: R^2_PGS <= h^2 <= {H_squared}")

    # Step 4: Evaluate each answer choice.
    print("\n--- Evaluating the Answer Choices ---")
    print(f"Based on the inequality 'R^2_PGS <= {H_squared}':")

    print("\nA. 'The polygenic score can not explain more than 50% of the variance...'")
    print("   This is NECESSARILY TRUE. The inequality R^2_PGS <= 0.5 directly supports this.")

    print("\nB. '...the PGS will approach a variance explained of 50%.'")
    print("   This is NOT NECESSARILY TRUE. This would only happen if h^2 = 0.5. It's possible that h^2 < 0.5 due to non-additive genetic effects, in which case the PGS would approach a value less than 0.5.")

    print("\nC. '...the PGS will not approach a variance explained of 50%...'")
    print("   This is NOT NECESSARILY TRUE. While often the case, it's theoretically possible that all genetic variance is additive (h^2 = 0.5). In that specific scenario, the PGS would approach 50%. The statement claims it 'will not', which is too strong.")

    print("\nD. '...epigenetic effects would limit the narrow-sense heritability to less than 0.5.'")
    print("   This is NOT NECESSARILY TRUE. The given H^2 of 0.5 already accounts for all sources of variance. The existence of epigenetic effects (which contribute to total variance) does not by itself necessitate that h^2 must be strictly less than 0.5.")

    print("\n--- Final Conclusion ---")
    print("The only statement that must be true based on the provided information is A.")

analyze_heritability_question()
<<<A>>>