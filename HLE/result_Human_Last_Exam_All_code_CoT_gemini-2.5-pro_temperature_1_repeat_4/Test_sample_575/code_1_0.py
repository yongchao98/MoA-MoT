def solve_heritability_problem():
    """
    This function explains the reasoning to solve the genetic heritability problem
    by demonstrating the mathematical relationships between variance components.
    """

    # --- Given Information ---
    # Vp: Total Phenotypic Variance
    # Vg: Total Genetic Variance
    # Ve: Total Environmental Variance
    # H2: Broad-sense Heritability
    #
    # The fundamental equation of quantitative genetics is Vp = Vg + Ve.
    # Broad-sense heritability (H2) is the proportion of phenotypic variance
    # explained by genetic variance.
    H2 = 0.5
    print(f"Given: Broad-sense heritability (H2 = Vg / Vp) is {H2}")
    print("-" * 30)

    # --- Analysis of the Polygenic Score (PGS) ---
    # A Polygenic Score (PGS) is a predictor of a phenotype based on genetic data.
    # Therefore, the variance it can explain (V_pgs) must be a component of,
    # and cannot exceed, the total genetic variance (Vg).
    print("Step 1: The variance explained by a polygenic score (V_pgs) is a subset of the total genetic variance (Vg).")
    print("This gives us the inequality: V_pgs <= Vg")
    print("-" * 30)

    # --- Relating PGS to Phenotypic Variance ---
    # We want to find the maximum proportion of the total phenotypic variance (Vp)
    # that the PGS can explain. This is the ratio V_pgs / Vp.
    # We can divide our inequality from Step 1 by Vp (which is a positive value).
    print("Step 2: To find the proportion of phenotypic variance explained by the PGS, we divide the inequality by Vp.")
    print("This gives us: (V_pgs / Vp) <= (Vg / Vp)")
    print("-" * 30)

    # --- Final Conclusion ---
    # We know from the problem statement that Vg / Vp is the broad-sense heritability, H2.
    print("Step 3: Substitute the given value of broad-sense heritability (H2) into the inequality.")
    print(f"We know that (Vg / Vp) = H2 = {H2}")
    print("\nTherefore, the final equation is:")
    
    # Final equation and conclusion
    # The variable names are kept for clarity in the output string
    proportion_explained_by_pgs = "V_pgs / Vp"
    print(f"{proportion_explained_by_pgs} <= {H2}")

    print("\nConclusion: This means the polygenic score can not explain more than 50% of the variance in the phenotype.")
    print("This directly corresponds to Answer Choice A.")

# Run the explanation
solve_heritability_problem()