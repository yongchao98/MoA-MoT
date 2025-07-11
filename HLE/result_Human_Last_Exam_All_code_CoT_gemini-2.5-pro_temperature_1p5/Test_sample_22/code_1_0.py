def explain_pgs_heritability():
    """
    Explains the relationship between SNP heritability and Polygenic Score (PGS) predictive ability.
    """

    # SNP heritability (h²_SNP) is the total proportion of phenotypic variance
    # explained by all considered SNPs. It is the theoretical maximum predictive power.
    explanation_h2 = "SNP heritability represents the theoretical ceiling of prediction."

    # A Polygenic Score (PGS) is a practical model built using estimated SNP effects
    # from a finite GWAS discovery sample. These estimates are not perfect.
    explanation_pgs = "A Polygenic Score is a practical tool built with imperfect, estimated data."

    # The variance explained by a PGS (R²) will be a fraction of the total
    # SNP heritability because the effect size estimates are noisy.
    # Therefore, R² from a PGS will be lower than h²_SNP.
    conclusion = "As a result, the predictive power of the PGS is necessarily lower than the total SNP heritability."

    final_answer = "True"

    print("Explanation:")
    print(f"1. {explanation_h2}")
    print(f"2. {explanation_pgs}")
    print(f"3. {conclusion}")
    print("\n------------------------------")
    print(f"The statement is: {final_answer}")

explain_pgs_heritability()