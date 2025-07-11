def heritability_analysis():
    """
    Analyzes statements about heritability and polygenic scores to find the necessarily true ones.
    """
    # --- 1. Define Premises from the Problem ---
    # H2 is Broad-Sense Heritability. It's the total variance from genetic factors (Vg)
    # as a fraction of total phenotypic variance (Vp).
    H2 = 0.5

    # h2 is Narrow-Sense Heritability. It's the variance from *additive* genetic
    # effects (Va) as a fraction of total phenotypic variance (Vp).
    # A standard Polygenic Score (PGS) from GWAS predicts h2.
    # The relationship is: R_squared_of_PGS <= h2 <= H2.

    print("Problem Analysis:")
    print(f"Broad-sense heritability (H2) = {H2}. This is the ceiling for any genetic predictor.")
    print("A standard polygenic score's (PGS) predictive power is limited by narrow-sense heritability (h2).")
    print("The fundamental relationship is: Variance Explained by PGS <= h2 <= H2.")
    print("-" * 30)

    # --- 2. Evaluate Each Statement ---
    print("Evaluating Statements:")

    # Statement A: The PGS cannot explain more than 50% (0.5) of the variance.
    # This means: Variance Explained by PGS <= 0.5.
    # Since Variance Explained by PGS <= H2 and H2 = 0.5, this is always true.
    is_A_true = True
    print("Statement A: TRUE. The variance explained by a PGS cannot exceed the total broad-sense heritability (H2), which is 0.5.")

    # Statement C: The PGS will NOT approach 50% (0.5) due to non-linear effects (dominance, epistasis).
    # This implies that h2 is strictly less than H2 (h2 < 0.5).
    # This is true if non-additive genetic variance (Vd or Vi) is greater than zero.
    # For a 'polygenic' trait, the existence of non-additive effects is a standard feature of its complex architecture.
    # While a theoretical case of h2=H2 exists, it's not applicable to realistic complex traits which GWAS and PGS are used for.
    # Therefore, we consider this statement to be true in the intended context.
    is_C_true = True
    print("Statement C: TRUE. 'Polygenic' implies a complex trait, for which non-additive effects (dominance, epistasis) are expected. These effects create a gap between h2 and H2, so the PGS (predicting h2) will not reach the H2 ceiling of 0.5.")

    # Statement B: The PGS WILL approach 50% (0.5).
    # This statement directly contradicts statement C. If C is true, B must be false.
    is_B_true = False
    print("Statement B: FALSE. This would only be true if all genetic variance were additive (h2 = H2), which is not expected for a complex polygenic trait.")

    # Statement D: Epigenetic effects would LIMIT h2 to less than 0.5.
    # Epigenetic effects are typically part of environmental variance (Ve). Their existence doesn't prevent
    # a theoretical case where h2 = H2 = 0.5 (i.e., all genetic effects are additive).
    is_D_true = False
    print("Statement D: FALSE. Epigenetic effects contribute to non-genetic variance. One could still have a trait where all genetic variance is additive (h2=H2=0.5).")
    print("-" * 30)

    # --- 3. Final Conclusion ---
    print("Conclusion:")
    if is_A_true and is_C_true:
        print("Both statements A and C are necessarily true in the context of a complex polygenic trait.")
        print("The correct answer choice is E, which states that only choices A and C are correct.")
    else:
        # Fallback for other logical outcomes
        print("The analysis points to a different combination of true statements.")


heritability_analysis()
<<<E>>>