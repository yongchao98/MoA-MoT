def analyze_heritability_statement():
    """
    Analyzes statements about a polygenic score based on given heritability.
    """

    # --- Step 1: Define known parameters and concepts from the problem ---
    # The problem states the broad-sense heritability (H^2) is 0.5.
    # H^2 represents the proportion of total phenotypic variance (Vp) that is
    # caused by all genetic variance (Vg).
    broad_sense_heritability = 0.5

    # Equation 1: H^2 = Vg / Vp
    print(f"The total variance explained by all genetic factors (broad-sense heritability) is {broad_sense_heritability}.")

    # Genetic variance (Vg) can be broken down into additive (Va), dominance (Vd),
    # and epistatic (Vi) components.
    # Equation 2: Vg = Va + Vd + Vi
    
    # Narrow-sense heritability (h^2) is the proportion of phenotypic variance
    # caused by only the additive genetic variance (Va).
    # Equation 3: h^2 = Va / Vp

    # A standard Polygenic Score (PGS) is an additive model built from GWAS effect sizes.
    # The maximum variance it can explain, even with perfect data, is the narrow-sense heritability (h^2).
    # Max variance explained by standard PGS = h^2

    print("From definitions, we know that additive variance (Va) cannot be greater than total genetic variance (Vg).")
    print("Therefore, narrow-sense heritability (h^2) must be less than or equal to broad-sense heritability (H^2).")
    print("-" * 60)

    # --- Step 2: Evaluate each statement ---

    # Statement A: The polygenic score can not explain more than 50% of the variance in the phenotype
    print("Analysis of Statement A:")
    print("A PGS is a predictor based on genetics. The absolute maximum variance any genetic predictor can explain is the total variance attributable to genetics, which is the broad-sense heritability (H^2).")
    print(f"Therefore, the variance explained by a PGS must be <= H^2.")
    # The equation representing this is: Variance explained by PGS <= h^2 <= H^2 = 0.5
    print(f"So, the PGS cannot explain more than {broad_sense_heritability*100}% of the variance. This statement is necessarily true.")
    print("\nConclusion for A: TRUE.\n")

    # Statement B: Given an arbitrarily large GWAS, the PGS will approach a variance explained of 50%
    print("Analysis of Statement B:")
    print("A perfect PGS (from an infinitely large GWAS) explains h^2 of the variance.")
    print(f"This would only be {broad_sense_heritability*100}% if h^2 = H^2, which means there are no non-additive effects (dominance or epistasis).")
    print("The problem does not state this, so h^2 could be less than 0.5. This statement is not necessarily true.")
    print("\nConclusion for B: FALSE.\n")

    # Statement C: Given an arbitrarily large GWAS, the PGS will not approach a variance explained of 50% due to non-linear effects
    print("Analysis of Statement C:")
    print("This statement claims that non-linear effects must exist, making h^2 strictly less than H^2.")
    print("While this is common for complex traits, the problem describes a 'theoretically ideal' population, where it is possible that all genetic variance is purely additive (h^2 = H^2).")
    print("Because it is not guaranteed that h^2 < H^2, this statement is not necessarily true.")
    print("\nConclusion for C: FALSE.\n")

    # Statement D: The existence of any epigenetic effects would limit the narrow-sense heritability to less than 0.5
    print("Analysis of Statement D:")
    print("Non-heritable epigenetic effects are typically categorized under environmental variance (Ve).")
    print("The existence of environmental variance does not require the existence of non-additive genetic variance (Vd or Vi).")
    print("A trait could have h^2 = H^2 = 0.5 even if environmental (including epigenetic) factors also contribute to the phenotype. This statement is not necessarily true.")
    print("\nConclusion for D: FALSE.\n")

    # --- Step 3: Conclude and output the final answer ---
    print("-" * 60)
    print("Based on the analysis, only Statement A is a necessary truth derived from the problem's premises.")
    final_answer = 'A'
    print(f"<<<{final_answer}>>>")

if __name__ == "__main__":
    analyze_heritability_statement()