def explain_heritability():
    """
    This script explains the relationship between broad-sense heritability,
    narrow-sense heritability, and the variance explained by a Polygenic Score (PGS).
    """

    print("--- Problem Analysis ---")
    print("We are given the following information:")
    H_sq_given = 0.5
    print(f"1. Broad-sense heritability (H²) = {H_sq_given}")
    print("2. A Polygenic Score (PGS) is built from typical GWAS data.")
    print("\n")

    print("--- Key Concepts in Quantitative Genetics ---")
    print("Vp = Total Phenotypic Variance")
    print("Vg = Total Genetic Variance")
    print("Ve = Environmental Variance")
    print("\nThe equation for Phenotypic Variance is: Vp = Vg + Ve")
    print("\nVg can be broken down into:")
    print("  Va = Additive Genetic Variance (effects of alleles sum up)")
    print("  Vd = Dominance Variance (interaction of alleles at the same locus)")
    print("  Vi = Epistatic Variance (interaction of alleles at different loci)")
    print("The equation for Genetic Variance is: Vg = Va + Vd + Vi")
    print("\n")

    print("--- Defining Heritability ---")
    print("Broad-sense heritability (H²) is the proportion of Vp due to ALL genetic factors.")
    print("H² = Vg / Vp = (Va + Vd + Vi) / Vp")
    print("\nNarrow-sense heritability (h²) is the proportion of Vp due to ADDITIVE genetic factors only.")
    print("h² = Va / Vp")
    print("\n")

    print("--- Connecting to the Polygenic Score (PGS) ---")
    print("A standard GWAS primarily detects additive effects. Therefore, a PGS built from it mainly captures Va.")
    print("The theoretical maximum variance a PGS can explain is equal to the narrow-sense heritability, h².")
    print("So, R²_pgs <= h²")
    print("\n")

    print("--- Evaluating the Statements ---")
    print(f"We know H² = Vg / Vp = {H_sq_given}")
    print("We also know Va <= Vg (since Vd and Vi cannot be negative).")
    print("Dividing by Vp gives: (Va / Vp) <= (Vg / Vp)")
    print("This means: h² <= H²")
    print(f"Since H² = {H_sq_given}, it must be true that h² <= {H_sq_given}.")
    print("\nConclusion:")
    print("The maximum variance a PGS can explain is h², and we've proven h² <= 0.5.")
    print("Therefore, the variance explained by the PGS cannot be more than 0.5 (or 50%).")
    print("\n")
    print("Let's re-examine the choices:")
    print("A: The PGS cannot explain more than 50% of the variance. This is NECESSARILY TRUE, as shown above.")
    print("B: The PGS will approach 50%. Not necessarily. This only happens if all genetic variance is additive (h² = H²).")
    print("C: The PGS will not approach 50%. Not necessarily. This is false if all genetic variance happens to be additive.")
    print("D: Epigenetics limits h² to < 0.5. Not necessarily. Heritable epigenetic effects could be purely additive, making h² = 0.5 possible.")
    print("\n")

    final_answer = "A"
    print(f"The only statement that is necessarily true is: {final_answer}")

explain_heritability()
<<<A>>>