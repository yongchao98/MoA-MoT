def analyze_heritability():
    """
    This function analyzes the relationship between heritability and Polygenic Scores (PGS)
    based on the premises given in the problem.
    """

    # --- GIVEN INFORMATION ---
    # Broad-sense heritability (H²) is the proportion of total phenotypic variance (Vp)
    # that is due to total genetic variance (Vg).
    H_squared = 0.5

    # --- DEFINITIONS ---
    # Vp = Vg + Ve (Phenotypic variance = Genetic variance + Environmental variance)
    # H² = Vg / Vp
    #
    # Vg = Va + Vd + Vi (Genetic variance = Additive + Dominance + Epistatic variance)
    #
    # Narrow-sense heritability (h²) = Va / Vp
    #
    # The variance explained by a Polygenic Score (R²_PGS) from a standard GWAS
    # is limited by the narrow-sense heritability (h²), which in turn is limited
    # by the broad-sense heritability (H²).

    # --- LOGICAL ANALYSIS ---

    print("Step 1: The problem states that broad-sense heritability (H²) is 0.5.")
    print("This means the total variance explained by all genetic factors (Vg) is 50% of the phenotypic variance (Vp).")
    print("Equation: Vg / Vp = 0.5")
    print("")

    print("Step 2: A Polygenic Score (PGS) is a model that uses only genetic information to predict a phenotype.")
    print("Therefore, the maximum variance a PGS can possibly explain (R²_PGS) is capped by the total genetic variance.")
    print("")

    print("Step 3: This establishes an inequality. The variance explained by the PGS must be less than or equal to the broad-sense heritability.")
    print("This is because R²_PGS is at most h² (Va/Vp), and h² is at most H² (Vg/Vp).")
    print("So, R²_PGS <= h² <= H².")
    print("")

    print("--- FINAL CONCLUSION ---")
    print("Based on the logic, we can form the final equation using the value provided:")
    # The final equation demonstrates the limit on the PGS's explanatory power.
    print(f"R²_PGS <= H² = {H_squared}")
    print("")
    print(f"This means the polygenic score can not explain more than {H_squared * 100}% of the variance in the phenotype.")
    print("This conclusion is robust and does not depend on the specific composition of the genetic variance (e.g., how much is additive vs. non-additive).")
    print("Therefore, statement A is the only choice that is necessarily true.")

analyze_heritability()
<<<A>>>