def analyze_heritability_statement():
    """
    Analyzes the provided statements about heritability and polygenic scores.
    This function will print the reasoning for each choice based on genetic principles.
    """

    # --- Problem Setup ---
    H_squared = 0.5  # Broad-sense heritability (Vg / Vp)
    phenotypic_variance_Vp = 100 # Assume a total variance for illustration

    # --- Derived Quantities ---
    # Vg is the total variance from ALL genetic factors.
    genetic_variance_Vg = H_squared * phenotypic_variance_Vp

    # The variance explained by a Polygenic Score (R_squared_PGS) is fundamentally
    # limited by the total contribution of genetics to the phenotype's variance.
    # The absolute ceiling for any genetic predictor is Vg.
    # Therefore, the proportion of variance explained by a PGS cannot exceed Vg / Vp.
    max_variance_explained_by_genetics = genetic_variance_Vg / phenotypic_variance_Vp

    print("--- Core Principles ---")
    print(f"Broad-sense heritability (H²) = {H_squared}")
    print(f"This means the total genetic variance (Vg) accounts for {H_squared*100}% of the total phenotypic variance (Vp).")
    print(f"A Polygenic Score (PGS) is a predictor based on genetic data.")
    print(f"Therefore, the maximum variance a PGS can possibly explain is limited by the total genetic variance.")
    print(f"Variance explained by PGS (R²_PGS) <= H²")
    print(f"R²_PGS <= {max_variance_explained_by_genetics}\n")


    print("--- Evaluating the Statements ---\n")

    # --- Statement A ---
    print("A. The polygenic score can not explain more than 50% of the variance in the phenotype.")
    print(f"   Reasoning: Based on our core principle, R²_PGS must be less than or equal to H², which is {H_squared}.")
    print("   This means the PGS can explain AT MOST 50% of the variance. It cannot explain more.")
    print("   Conclusion: This statement is NECESSARILY TRUE.\n")

    # --- Statement B & C ---
    # Narrow-sense heritability (h²) = Va / Vp.
    # A standard, linear PGS approaches h² as GWAS size increases.
    # We only know that h² <= H², because Va <= Vg. We don't know if they are equal.
    print("B. Given an arbitrarily large GWAS, the polygenic score will approach a variance explained of 50%.")
    print("C. Given an arbitrarily large GWAS, the polygenic score ... will not approach ... 50% due to non-linear effects.")
    print("   Reasoning for B & C: A standard PGS approaches narrow-sense heritability (h²), not broad-sense (H²).")
    print("   We know h² <= H² (0.5).")
    print("   - If all genetic variance is additive (no dominance/epistasis), then h² = H² = 0.5, making B true and C false.")
    print("   - If there is any non-additive variance, then h² < H², making C true and B false.")
    print("   Since we don't know the composition of the genetic variance, neither B nor C is *necessarily* true.\n")

    # --- Statement D ---
    print("D. The existence of any epigenetic effects would limit the narrow-sense heritability to less than 0.5.")
    print("   Reasoning: Epigenetic effects are typically considered part of the environmental variance (Ve).")
    print("   The composition of Ve does not constrain the composition of Vg. It's possible for h² to equal H² (0.5) even if epigenetic effects exist.")
    print("   Conclusion: This statement is NOT necessarily true.\n")

    # --- Statement E ---
    print("E. None of the other answer choices are correct.")
    print("   Reasoning: Statement A is necessarily true.")
    print("   Conclusion: This statement is FALSE.\n")


# Execute the analysis
analyze_heritability_statement()
