def analyze_heritability_and_pgs():
    """
    Analyzes the relationship between heritability and the predictive power of a polygenic score (PGS).
    This script will demonstrate why only statement A is necessarily true.
    """
    # Step 1: Define the core concepts and the given value from the problem.
    H2_broad_sense = 0.5

    print("--- Core Concepts & Definitions ---")
    print(f"Broad-sense heritability (H²), the proportion of phenotypic variance (Vp) due to ALL genetic variance (Vg), is given as: {H2_broad_sense}")
    print("H² = Vg / Vp")
    print("\nGenetic variance (Vg) can be decomposed into: Vg = Va + Vd + Vi")
    print("  - Va: Additive genetic variance (captured by a standard PGS)")
    print("  - Vd & Vi: Non-additive variance (dominance and epistasis)")
    print("\nNarrow-sense heritability (h²) is the proportion of phenotypic variance due to ADDITIVE genetic variance:")
    print("h² = Va / Vp")
    print("\nThe maximum variance a standard PGS can explain is h².")
    print("-" * 45)

    # Step 2: Perform the logical deduction.
    print("\n--- Logical Deduction ---")
    print("We are given the equation: Vg / Vp = 0.5")
    print("By definition, we know that Va must be less than or equal to Vg, so: Va <= Vg")
    print("Dividing the inequality by Vp gives: Va / Vp <= Vg / Vp")
    print("Substituting the heritability definitions, we get: h² <= H²")
    print(f"Therefore, h² must be less than or equal to {H2_broad_sense}.")
    print("\nConclusion: The variance explained by a PGS (at most h²) cannot be more than 0.5 (or 50%).")
    print("-" * 45)

    # Step 3: Use scenarios to test the answer choices.
    # Let's assume a total phenotypic variance (Vp) of 100 for easy calculation.
    Vp = 100.0
    Vg = H2_broad_sense * Vp
    print(f"\n--- Demonstrating with Scenarios (assuming Vp = {Vp}) ---")
    print(f"Given H² = {H2_broad_sense}, the total genetic variance Vg must be {Vg}.")

    # Scenario 1: All genetic variance is additive. This tests statement C.
    print("\n>>> Scenario 1: All genetic variance is additive (no dominance or epistasis).")
    Va_1 = 50.0
    Vd_1, Vi_1 = 0.0, 0.0
    h2_1 = Va_1 / Vp
    print(f"Let's assume: Va={Va_1}, Vd={Vd_1}, Vi={Vi_1}. Total Vg = {Va_1 + Vd_1 + Vi_1}")
    print(f"The resulting equation for h² is: Va / Vp = {Va_1} / {Vp} = {h2_1}")
    print("In this case, the PGS explained variance (h²) is 50%. This shows that statement C (that it CAN'T be 50%) is not necessarily true.")

    # Scenario 2: Some genetic variance is non-additive. This tests statement B.
    print("\n>>> Scenario 2: Genetic variance includes non-additive components.")
    Va_2 = 40.0
    Vd_2, Vi_2 = 5.0, 5.0
    h2_2 = Va_2 / Vp
    print(f"Let's assume: Va={Va_2}, Vd={Vd_2}, Vi={Vi_2}. Total Vg = {Va_2 + Vd_2 + Vi_2}")
    print(f"The resulting equation for h² is: Va / Vp = {Va_2} / {Vp} = {h2_2}")
    print("In this case, the PGS explained variance (h²) is 40%. This shows that statement B (that it WILL be 50%) is not necessarily true.")

    # Final conclusion based on the logic that holds in all scenarios.
    print("\n" + "-" * 45)
    print("\nFinal Conclusion:")
    print("In all possible scenarios, the PGS variance explained (h²) is less than or equal to H² (0.5).")
    print("Therefore, the only statement that is NECESSARILY true is A.")

analyze_heritability_and_pgs()
