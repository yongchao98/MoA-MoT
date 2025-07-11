def analyze_heritability_problem():
    """
    This script explains the reasoning behind the correct answer choice
    by defining the key concepts and evaluating each option.
    """

    # --- Definitions and Given Information ---
    H2 = 0.5  # Broad-sense heritability (Vg / Vp)

    print("Problem Analysis:")
    print("=" * 25)
    print(f"Given: Broad-sense heritability (H²) = {H2}")
    print("This means the total variance explained by all genetic factors is 50%.")
    print("The remaining 50% is due to environmental factors.\n")

    print("Key Concepts:")
    print("- A Polygenic Score (PGS) is a predictor based only on genetic data.")
    print("- The maximum possible variance any genetic predictor can explain is H².\n")

    # --- Evaluation of Answer Choices ---
    print("Evaluation:")
    print("-" * 25)

    # Statement A
    print("A: 'The PGS cannot explain more than 50% of the variance.'")
    print(f"   - Reasoning: The performance of any genetic predictor is capped by the total genetic variance, H² ({H2*100}%).")
    print("   - Verdict: This is necessarily TRUE.\n")

    # Statement B
    print("B: 'The PGS will approach a variance explained of 50%.'")
    print("   - Reasoning: A standard PGS captures additive genetic variance (h²), not total genetic variance (H²).")
    print("   - Since h² <= H², the PGS will approach a value <= 50%, not necessarily 50%.")
    print("   - Verdict: Not necessarily true.\n")

    # Statement C
    print("C: 'The PGS will not approach 50% due to non-linear effects.'")
    print("   - Reasoning: This assumes non-linear effects (dominance, epistasis) must exist.")
    print("   - In a 'theoretically ideal' case, all genetic variance could be additive (h² = H² = 0.5).")
    print("   - In that case, the PGS would approach 50%.")
    print("   - Verdict: Not necessarily true.\n")

    # Statement D
    print("D: 'Epigenetic effects would limit narrow-sense heritability to less than 0.5.'")
    print("   - Reasoning: Epigenetic effects contribute to genetic variance but don't require non-additive variance.")
    print("   - It's still theoretically possible for h² = H² = 0.5 even with epigenetic effects.")
    print("   - Verdict: Not necessarily true.\n")

    # --- Final Conclusion ---
    final_answer = 'A'
    print("Conclusion:")
    print("=" * 25)
    print("Statement A is the only choice that is a direct, definitional consequence of H²=0.5.")
    
# Execute the analysis
analyze_heritability_problem()