def analyze_heritability_statements():
    """
    Analyzes statements about heritability and polygenic scores based on given parameters.
    This function prints a step-by-step evaluation of each statement.
    """
    # --- Given Information ---
    H2 = 0.5  # Broad-sense heritability

    print("--- Foundational Equations ---")
    print("Phenotypic Variance (Vp) = Genetic Variance (Vg) + Environmental Variance (Ve)")
    print("Genetic Variance (Vg) = Additive (Va) + Dominance (Vd) + Epistatic (Vi)")
    print(f"Broad-Sense Heritability (H2) = Vg / Vp = {H2}")
    print("Narrow-Sense Heritability (h2) = Va / Vp")
    print("Variance explained by an ideal Polygenic Score (R2_PGS) approaches h2.")
    print("-" * 40)

    # --- Statement Evaluation ---
    print("\n--- Evaluating Answer Choices ---")

    # A. The polygenic score can not explain more than 50% of the variance in the phenotype
    print("\n[A] Analysis:")
    print("R2_PGS is limited by h2.")
    print("h2 = Va / Vp.")
    print(f"H2 = (Va + Vd + Vi) / Vp = {H2}.")
    print("Since Va <= (Va + Vd + Vi), it follows that h2 <= H2.")
    print(f"Therefore, R2_PGS <= h2 <= {H2}.")
    print("Conclusion: Statement A is TRUE. The PGS cannot explain more than 50% of the variance.")

    # B. Given an arbitrarily large GWAS, the polygenic score will approach a variance explained of 50%
    print("\n[B] Analysis:")
    print("This implies that R2_PGS will approach H2, meaning h2 = H2.")
    print(f"For h2 to equal H2 ({H2}), the non-additive genetic variance (Vd + Vi) must be 0.")
    print("This is the special case of a purely additive trait, which is not guaranteed.")
    print("Conclusion: Statement B is NOT necessarily true.")

    # C. Given an arbitrarily large GWAS, the polygenic score ... will not approach ... 50% due to ... non-linear effects
    print("\n[C] Analysis:")
    print("This statement claims R2_PGS will not approach 50% because of non-additive effects (dominance and epistasis).")
    print("This means h2 < H2 because Vd + Vi > 0.")
    print("For a 'polygenic' trait, it is a standard biological assumption that non-additive effects exist, so Vd + Vi > 0.")
    print(f"This makes h2 < H2 ({H2}).")
    print("Conclusion: Statement C is TRUE under standard biological assumptions.")

    # D. The existence of any epigenetic effects would limit the narrow-sense heritability to less than 0.5
    print("\n[D] Analysis:")
    print("Non-heritable epigenetic effects contribute to Environmental Variance (Ve).")
    print("The value of H2 = 0.5 already accounts for the total Vp = Vg + Ve.")
    print("It is possible for a trait to be purely additive (h2 = H2 = 0.5) while still having Ve > 0.")
    print("The existence of Ve does not force h2 to be strictly less than H2.")
    print("Conclusion: Statement D is NOT necessarily true.")

    print("-" * 40)
    print("\nFinal determination: Statements A and C are the correct assertions.")

# Run the analysis
analyze_heritability_statements()