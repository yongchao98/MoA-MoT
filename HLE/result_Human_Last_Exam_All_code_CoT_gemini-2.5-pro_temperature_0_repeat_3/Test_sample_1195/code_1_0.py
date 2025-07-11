def solve_genetics_cross():
    """
    Calculates and explains the F2 phenotypic ratio for a Drosophila cross
    involving an X-linked gene (vermilion) and an autosomal suppressor gene (su-v).
    """
    # Step 1: Define Parental (P) Genotypes and Phenotypes
    # v = vermilion allele (X-linked, recessive)
    # su-v = suppressor allele (autosomal, recessive)
    # su-v+ = wild-type suppressor allele (dominant)
    # Interaction: su-v/su-v genotype restores wild-type eyes in flies with the vermilion genotype.
    p_female_genotype = "X(v)X(v) su-v/su-v"
    p_male_genotype = "X(v)Y su-v+/su-v+"

    print("--- Step 1: Parental (P) Generation ---")
    print(f"Female Genotype: {p_female_genotype}")
    print(f"Male Genotype:   {p_male_genotype}")
    print("\n")

    # Step 2: Determine F1 Generation
    # F1 Cross: X(v)X(v) su-v/su-v  x  X(v)Y su-v+/su-v+
    # All F1 offspring are heterozygous for the suppressor gene (su-v+/su-v)
    # and have the vermilion allele (v). Since su-v is recessive, all F1s are vermilion.
    print("--- Step 2: F1 Generation ---")
    print("All F1 offspring have the genotype X(v)X(v) su-v+/su-v (female) or X(v)Y su-v+/su-v (male).")
    print("Because the suppressor allele (su-v) is recessive, all F1 offspring have vermilion eyes.")
    print("\n")

    # Step 3: Analyze F2 Generation from F1 Intercross
    # F1 Cross: X(v)X(v) su-v+/su-v  x  X(v)Y su-v+/su-v
    print("--- Step 3: F2 Generation Analysis ---")
    print("The F2 generation results from intercrossing the F1 flies.")
    print("1. Analysis of the X-linked gene (v):")
    print("  The cross is X(v)X(v) x X(v)Y. All F2 offspring will have the genetic basis for vermilion eyes.")
    print("\n2. Analysis of the autosomal suppressor gene (su-v):")
    print("  The cross is su-v+/su-v x su-v+/su-v. This gives a genotypic ratio of:")
    print("  1/4 su-v+/su-v+ : 2/4 su-v+/su-v : 1/4 su-v/su-v")
    print("\n")

    # Step 4: Combine results to find the F2 phenotypic ratio
    # The phenotype is determined by the su-v genotype, since all F2s have the v allele.
    prob_no_suppression_numerator = 1 + 2
    prob_suppression_numerator = 1
    denominator = 4

    print("--- Step 4: Final F2 Phenotypic Ratio ---")
    print("Combining the results:")
    print(f" - Offspring with at least one su-v+ allele ({prob_no_suppression_numerator}/{denominator} of total) will be VERMILION.")
    print(f" - Offspring with the su-v/su-v genotype ({prob_suppression_numerator}/{denominator} of total) will be WILD-TYPE (suppressed).")
    print("\nTherefore, the final equation for the phenotypic ratio is:")
    print(f"{prob_no_suppression_numerator}/{denominator} vermilion : {prob_suppression_numerator}/{denominator} wild-type")

solve_genetics_cross()