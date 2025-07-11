def analyze_heritability_statement():
    """
    Analyzes statements about heritability and polygenic scores.
    """
    # --- Step 1 & 2: Define terms and relationships ---
    H2 = 0.5  # Given Broad-Sense Heritability

    print("--- Understanding the Concepts ---")
    print(f"Broad-Sense Heritability (H²): {H2}")
    print("H² represents the proportion of total phenotypic variance (Vp) that is due to ALL genetic factors (Vg).")
    print("Equation: H² = Vg / Vp")
    print("Genetic Variance (Vg) can be broken down: Vg = Va + Vd + Vi")
    print("  - Va: Additive genetic variance (effects of alleles sum up linearly)")
    print("  - Vd: Dominance variance (non-linear interactions between alleles at the same locus)")
    print("  - Vi: Epistatic variance (non-linear interactions between alleles at different loci)")
    print("\nNarrow-Sense Heritability (h²):")
    print("h² represents the proportion of phenotypic variance due to ADDITIVE genetic factors only.")
    print("Equation: h² = Va / Vp")
    print("\nPolygenic Score (PGS):")
    print("A standard PGS is built from GWAS summary statistics, which primarily capture additive effects.")
    print("Therefore, the maximum variance a standard PGS can explain is h².")
    print("\nKey Relationship: Since Vg = Va + Vd + Vi, it is always true that Va <= Vg.")
    print("This means that h² <= H². In our case, h² <= 0.5.")
    print("-" * 30)

    # --- Step 3 & 4: Analyze each answer choice ---
    print("\n--- Evaluating the Answer Choices ---")

    # Let's assume a total phenotypic variance (Vp) of 100 units for easy calculation.
    # Vp = 100
    # Vg = H² * Vp = 0.5 * 100 = 50
    # Vg = Va + Vd + Vi = 50

    print("\nA. The polygenic score can not explain more than 50% of the variance in the phenotype.")
    print("   - The total variance from ALL genetic sources is H² = 50%.")
    print("   - A PGS is a predictor based on genetic information.")
    print("   - By definition, no genetic predictor can explain more variance than the total variance caused by genetics.")
    print("   - Therefore, the PGS explained variance must be <= H² (0.5).")
    print("   - This statement is NECESSARILY TRUE.")

    print("\nB. Given an arbitrarily large GWAS, the polygenic score will approach a variance explained of 50%.")
    print("   - A PGS from a large GWAS will approach explaining a variance of h² (Va / Vp).")
    print("   - This statement claims h² will approach 0.5 (or H²).")
    print("   - This is only true if all genetic variance is additive (Vd = 0 and Vi = 0).")
    print("   - Counterexample: If Vg=50, but Va=40, Vd=5, Vi=5, then H²=0.5 but h²=0.4.")
    print("   - In this case, the PGS would approach 40% variance explained, not 50%.")
    print("   - This statement is NOT NECESSARILY TRUE.")

    print("\nC. Given an arbitrarily large GWAS, the polygenic score ... will not approach ... 50% due to non-linear effects...")
    print("   - This statement claims that non-linear effects (Vd or Vi) MUST exist, so h² < H².")
    print("   - But it's theoretically possible that for this trait, all genetic variance is additive.")
    print("   - Counterexample: If Vg=50, and it's all additive (Va=50, Vd=0, Vi=0), then H²=0.5 and h²=0.5.")
    print("   - In this case, the PGS would approach 50% variance explained.")
    print("   - This statement is NOT NECESSARILY TRUE.")

    print("\nD. The existence of any epigenetic effects would limit the narrow-sense heritability to less than 0.5.")
    print("   - Epigenetic effects are typically considered part of the environmental variance (Ve) or gene-environment interactions.")
    print("   - The broad-sense heritability (H²=0.5) is a given ratio. This value already accounts for the total phenotypic variance, whatever its cause.")
    print("   - The definitional relationship h² <= H² (so h² <= 0.5) always holds, regardless of what constitutes the environmental variance.")
    print("   - The statement is not necessarily true; it doesn't change the fundamental relationship.")

    print("\nE. None of the other answer choices are correct.")
    print("   - Based on the analysis above, statement A is necessarily true.")
    print("-" * 30)

    # --- Step 5: Final Conclusion ---
    print("\nConclusion: The only statement that is necessarily true is A.")
    print("The total genetic variance provides a hard ceiling (H² = 0.5) on the variance that can be explained by any genetic predictor, including a polygenic score.")

analyze_heritability_statement()
<<<A>>>