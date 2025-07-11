def analyze_heritability_problem():
    """
    Analyzes the relationship between heritability and polygenic scores
    based on the provided information.
    """

    # Given information from the problem
    broad_sense_heritability_H2 = 0.5

    print("--- Problem Analysis ---")
    print(f"Given: Broad-sense heritability (H²) = {broad_sense_heritability_H2}")
    print("\n--- Key Definitions ---")
    print("Vp = Total Phenotypic Variance")
    print("Vg = Total Genetic Variance (Va + Vd + Vi)")
    print("Va = Additive Genetic Variance")
    print("\nH² = Vg / Vp  (Broad-sense heritability)")
    print("h² = Va / Vp  (Narrow-sense heritability)")
    print("\n--- Core Logic ---")
    print("1. By definition, Total Genetic Variance (Vg) is the sum of Additive (Va), Dominance (Vd), and Epistatic (Vi) variance.")
    print("2. Since variance cannot be negative, Va must be less than or equal to Vg (Va <= Vg).")
    print("3. This implies that narrow-sense heritability (h²) must be less than or equal to broad-sense heritability (H²).")
    print("   Equation: h² <= H²")
    print("4. A standard Polygenic Score (PGS) is built from additive effects, so its maximum explanatory power is h².")
    print("   Equation: Maximum R² of PGS = h²")
    print("\n--- Evaluating Statement A ---")
    print("Statement A: The polygenic score can not explain more than 50% of the variance.")
    print("Let's test this: Is `Maximum R² of PGS <= 0.5` a necessary truth?")
    print(f"From step 4, we substitute h²:  h² <= {broad_sense_heritability_H2}")
    print(f"From step 3, we know h² <= H²: h² <= {broad_sense_heritability_H2}")
    print(f"The final equation is: Maximum R² of PGS <= {broad_sense_heritability_H2}")
    print("\nConclusion: The statement is necessarily true. The maximum variance a PGS can explain is h², which cannot be greater than H² (0.5).")

analyze_heritability_problem()