def analyze_heritability():
    """
    Analyzes the relationship between broad-sense heritability, narrow-sense heritability,
    and the predictive power of a polygenic score (PGS).
    """

    # --- Given Information ---
    # Broad-sense heritability (H_sq) is the proportion of total phenotypic variance (Vp)
    # that is due to all genetic variance (Vg).
    H_sq = 0.5
    print(f"Given: Broad-sense heritability (H²) = Vg / Vp = {H_sq}")

    # --- Definitions ---
    # Genetic variance (Vg) has several components:
    # Va: Additive genetic variance (effects of alleles summed up)
    # Vd: Dominance variance (interaction between alleles at the same locus)
    # Vi: Epistatic variance (interaction between alleles at different loci)
    print("Definition: Total Genetic Variance (Vg) = Va + Vd + Vi")

    # Narrow-sense heritability (h_sq) is the proportion of phenotypic variance
    # due to ONLY the additive genetic variance (Va).
    print("Definition: Narrow-sense heritability (h²) = Va / Vp")

    # A standard Polygenic Score (PGS) from GWAS is a linear sum of SNP effects,
    # so it primarily captures additive genetic variance.
    # Therefore, the maximum variance a PGS can explain is h².
    print("Fact: Max variance explained by a standard PGS is h².")

    # --- Logical Deduction ---
    print("\n--- Logical Deduction ---")

    # Since Vd and Vi are variances, they must be greater than or equal to zero.
    # Vd >= 0, Vi >= 0
    # Therefore, the additive variance (Va) must be less than or equal to the total genetic variance (Vg).
    print("Step 1: Because Vg = Va + Vd + Vi, and Vd >= 0, Vi >= 0, it must be that Va <= Vg.")

    # Divide the inequality by the total phenotypic variance (Vp).
    # The inequality relationship still holds.
    print("Step 2: Dividing by Vp, we get (Va / Vp) <= (Vg / Vp).")

    # Substitute the definitions of h² and H².
    # (Va / Vp) is h²
    # (Vg / Vp) is H²
    print("Step 3: By definition, this means h² <= H².")

    # Apply the given value for H².
    print(f"Step 4: Substituting the given H² = {H_sq}, we find that h² <= {H_sq}.")

    # --- Conclusion ---
    print("\n--- Conclusion ---")
    print(f"The maximum variance a polygenic score can explain is h².")
    print(f"We have proven that h² must be less than or equal to {H_sq} (50%).")
    print("Therefore, the polygenic score cannot explain MORE than 50% of the variance.")
    print("This means statement A is necessarily true.")

analyze_heritability()
<<<A>>>