def analyze_heritability_statements():
    """
    Analyzes statements about heritability and polygenic scores based on given information.
    """
    # Given: Broad-sense heritability (H²) is 0.5.
    H_squared = 0.5

    # This means the proportion of phenotypic variance (Vp) due to total genetic variance (Vg) is 0.5.
    # H² = Vg / Vp = 0.5

    # Total genetic variance (Vg) consists of additive (Va), dominance (Vd), and epistatic (Vi) components.
    # Vg = Va + Vd + Vi

    # A Polygenic Score (PGS) is a predictor based on genetic information.
    # The variance explained by any genetic predictor (PGS_R_squared) cannot exceed the total
    # variance contributed by all genetic factors (H²).
    # This establishes a fundamental upper limit.

    print("--- Evaluating the Statements ---")
    print(f"Given: Broad-sense heritability (H²) = {H_squared}\n")

    print("Statement A: The polygenic score can not explain more than 50% of the variance in the phenotype.")
    print("Reasoning:")
    print("The variance explained by a PGS (PGS_R²) is a part of the total genetic variance (Vg).")
    print("Therefore, the proportion of phenotypic variance explained by a PGS cannot exceed the broad-sense heritability (H²).")
    print("The key relationship is:")
    print(f"PGS_R² <= H²")
    print(f"Substituting the given value: PGS_R² <= {H_squared}")
    print("This statement is a direct consequence of the definition of broad-sense heritability. It is NECESSARILY TRUE.")
    print("-" * 40)

    print("Statement B: Given an arbitrarily large GWAS, the polygenic score will approach a variance explained of 50%.")
    print("Reasoning:")
    print("A standard PGS from GWAS captures additive effects, so its explained variance approaches narrow-sense heritability (h² = Va / Vp).")
    print("We know h² <= H². h² only equals H² if dominance and epistatic variances are zero.")
    print("Since we don't know if this is the case, we cannot be sure the PGS will approach 0.5. This statement is NOT necessarily true.")
    print("-" * 40)

    print("Statement C: ...the PGS...will not approach a variance explained of 50% due to gene-gene interactions and other non-linear effects...")
    print("Reasoning:")
    print("This statement assumes that non-linear effects (dominance, epistasis) MUST exist.")
    print("While common, it is theoretically possible for a trait to have only additive genetic variance. In that case, h² = H² = 0.5, and the PGS would approach 0.5.")
    print("Since we cannot rule out this possibility, this statement is NOT necessarily true.")
    print("-" * 40)
    
    print("Statement D: The existence of any epigenetic effects would limit the narrow-sense heritability to less than 0.5.")
    print("Reasoning:")
    print("By definition, h² <= H², so h² is already <= 0.5. Epigenetic effects (often considered part of environmental variance) don't cause this definitional inequality.")
    print("This statement is NOT necessarily true.")
    print("-" * 40)

if __name__ == '__main__':
    analyze_heritability_statements()
    print("\n<<<A>>>")