def analyze_heritability():
    """
    Illustrates the concepts of heritability to analyze the problem statements.

    This code demonstrates two scenarios based on the given information:
    1. A hypothetical case with only additive genetic variance.
    2. A realistic case with both additive and non-additive genetic variance.
    """
    # --- GIVENS AND DEFINITIONS ---
    # H2: Broad-sense heritability (Vg / Vp)
    # h2: Narrow-sense heritability (Va / Vp)
    # Vp: Total Phenotypic Variance
    # Vg: Total Genetic Variance = Va (Additive) + Vd (Dominance) + Vi (Epistatic)
    # R2_PGS: Variance explained by an ideal linear Polygenic Score, which equals h2.

    H2 = 0.5
    # For demonstration, we assume a total phenotypic variance of 100 units.
    Vp = 100.0
    # From the broad-sense heritability, we calculate the total genetic variance.
    Vg = H2 * Vp

    print("### Problem Analysis ###")
    print(f"Given Broad-Sense Heritability (H\u00b2) = {H2}")
    print(f"Let's assume a Total Phenotypic Variance (Vp) = {Vp}")
    print(f"This means the Total Genetic Variance (Vg) = H\u00b2 * Vp = {H2} * {Vp} = {Vg}\n")

    # --- SCENARIO 1: Purely Additive Genetics (Hypothetical) ---
    print("--- Scenario 1: Purely Additive Genetic Model ---")
    Va_1 = Vg  # All genetic variance is additive
    Vd_plus_Vi_1 = 0.0
    h2_1 = Va_1 / Vp
    R2_PGS_1 = h2_1 * 100

    print(f"In this case, Va = {Va_1}, and non-additive variance (Vd + Vi) = {Vd_plus_Vi_1}")
    print(f"The Narrow-Sense Heritability (h\u00b2) = Va / Vp = {Va_1} / {Vp} = {h2_1}")
    print(f"An ideal PGS would explain {R2_PGS_1:.1f}% of phenotypic variance.")
    print("Conclusion: Here, the PGS *does* approach 50%. This shows statement C is not true in *all* logical possibilities.\n")


    # --- SCENARIO 2: Genetics with Non-Additive Effects (Realistic) ---
    print("--- Scenario 2: Realistic Model with Non-Additive Effects ---")
    # Let's assume 70% of genetic variance is additive, and 30% is non-additive.
    Va_2 = Vg * 0.7
    Vd_plus_Vi_2 = Vg * 0.3
    h2_2 = Va_2 / Vp
    R2_PGS_2 = h2_2 * 100

    print(f"In this case, Va = {Va_2}, and non-additive variance (Vd + Vi) = {Vd_plus_Vi_2}")
    print(f"The Narrow-Sense Heritability (h\u00b2) = Va / Vp = {Va_2} / {Vp} = {h2_2}")
    print(f"An ideal PGS would explain {R2_PGS_2:.1f}% of phenotypic variance.")
    print("Conclusion: The PGS variance explained is less than 50% due to non-additive effects. This reflects the common understanding of polygenic traits and makes statement C true in context.\n")

    # --- FINAL EVALUATION ---
    print("### Overall Conclusion ###")
    print("Statement A: 'The PGS cannot explain more than 50% of the variance.'")
    print("This is ALWAYS TRUE, because the max variance explained is h\u00b2, and h\u00b2 \u2264 H\u00b2 (which is 0.5).\n")

    print("Statement C: 'The PGS will not approach 50% due to non-linear effects.'")
    print("This is TRUE in any realistic scenario for a complex polygenic trait (like Scenario 2), which is the context of the question.\n")
    
    print("Therefore, both A and C are considered correct statements.")

if __name__ == '__main__':
    analyze_heritability()