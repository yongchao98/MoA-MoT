def explain_heritability():
    """
    Explains the reasoning behind the correct answer based on quantitative genetics principles.
    """
    
    # The problem provides the broad-sense heritability (H^2).
    H_squared = 0.5

    # Broad-sense heritability (H^2) is the proportion of phenotypic variance (Vp)
    # explained by total genetic variance (Vg).
    # Equation: H^2 = Vg / Vp

    # Narrow-sense heritability (h^2) is the proportion of phenotypic variance (Vp)
    # explained by *additive* genetic variance (Va).
    # Equation: h^2 = Va / Vp

    # A polygenic score (PGS) based on a standard GWAS primarily captures additive genetic effects.
    # Therefore, the maximum variance a PGS can explain is equal to the narrow-sense heritability (h^2).
    # Let R2_PGS be the variance explained by the PGS.
    # So, R2_PGS_max = h^2

    # Total genetic variance (Vg) is the sum of additive (Va) and non-additive (V_non-additive) components.
    # Equation: Vg = Va + V_non-additive
    # This means that Va is always less than or equal to Vg (Va <= Vg).

    # From this relationship, we can derive the relationship between h^2 and H^2:
    # Va / Vp <= Vg / Vp
    # h^2 <= H^2

    # Now, we apply the number from the problem.
    # We are given that H^2 = 0.5.
    # The variance explained by the PGS (R2_PGS) is limited by h^2.
    # The final equation is: R2_PGS <= h^2 <= H^2
    
    print("The final relationship is: (Variance explained by PGS) <= (Narrow-sense heritability) <= (Broad-sense heritability)")
    print(f"Based on the numbers in the problem, this means: (Variance explained by PGS) <= {H_squared}")
    print("\nStatement A says the polygenic score cannot explain more than 50% of the variance.")
    print("This is equivalent to saying: (Variance explained by PGS) <= 0.5, which is necessarily true.")

explain_heritability()
<<<A>>>