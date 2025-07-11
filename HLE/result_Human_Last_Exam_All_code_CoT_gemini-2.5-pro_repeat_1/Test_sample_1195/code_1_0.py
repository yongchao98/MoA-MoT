def solve_genetics_cross():
    """
    Solves the Drosophila gene interaction problem step-by-step.
    """
    print("### Solving the Drosophila Gene Interaction Problem ###\n")

    # Step 1: Define Parental (P) Generation
    print("Step 1: Define Parental (P) Generation Genotypes")
    # Note: Xv represents the X chromosome with the vermilion allele. su-v+ is the wild-type, non-suppressing allele.
    p_female_genotype = "X(v)X(v); su-v/su-v"
    p_male_genotype = "X(v)Y; su-v+/su-v+"
    print(f"Parental Female Genotype: {p_female_genotype}")
    print(f"Parental Male Genotype: {p_male_genotype}\n")

    # Step 2: Determine F1 Generation
    print("Step 2: Determine F1 Generation")
    print(" - The P female produces only one type of gamete: X(v); su-v")
    print(" - The P male produces two types of gametes: X(v); su-v+ and Y; su-v+")
    f1_female_genotype = "X(v)X(v); su-v+/su-v"
    f1_male_genotype = "X(v)Y; su-v+/su-v"
    print(f"Resulting F1 Female Genotype: {f1_female_genotype}")
    print(f"Resulting F1 Male Genotype: {f1_male_genotype}")
    print("Phenotype: All F1 flies are heterozygous for the suppressor gene (su-v+/su-v), so no suppression occurs. Since all F1s also have the vermilion genotype, they all have vermilion eyes.\n")

    # Step 3: Analyze the F1 x F1 Cross to get the F2 Generation
    print("Step 3: Analyze the F1 x F1 Cross for the F2 Generation")
    print("We can analyze the two genes independently.\n")
    
    print("  A) Autosomal Cross (su-v+/su-v x su-v+/su-v):")
    no_suppression_num = 3
    no_suppression_den = 4
    suppression_num = 1
    suppression_den = 4
    print(f"     - Genotypic Ratio: 1/4 su-v+/su-v+ : 2/4 su-v+/su-v : 1/4 su-v/su-v")
    print(f"     - This means {no_suppression_num}/{no_suppression_den} of F2 offspring have a non-suppressing genotype (su-v+/_).")
    print(f"     - And {suppression_num}/{suppression_den} of F2 offspring have the recessive suppressing genotype (su-v/su-v).\n")

    print("  B) X-linked Cross (X(v)X(v) x X(v)Y):")
    print("     - All offspring receive an X(v) from the mother.")
    print("     - Therefore, 100% of the F2 generation has the vermilion genotype (X(v)X(v) or X(v)Y).\n")

    # Step 4: Calculate Final F2 Phenotypic Ratio
    print("Step 4: Calculate the Final F2 Phenotypic Ratio")
    print("The phenotype is determined by whether the vermilion genotype is suppressed or not.\n")

    print("  - Vermilion Phenotype:")
    print("    Condition: Vermilion genotype (100%) AND Non-suppressing genotype (su-v+/_).")
    print(f"    Calculation: 1 * ({no_suppression_num}/{no_suppression_den}) = {no_suppression_num}/{no_suppression_den}")
    
    print("\n  - Wild-type Phenotype:")
    print("    Condition: Vermilion genotype (100%) AND Suppressing genotype (su-v/su-v).")
    print(f"    Calculation: 1 * ({suppression_num}/{suppression_den}) = {suppression_num}/{suppression_den}")

    print("\n-------------------------------------------------")
    print("Final Expected F2 Phenotypic Ratio:")
    print(f"{no_suppression_num}/{no_suppression_den} vermilion : {suppression_num}/{suppression_den} wild-type")
    print("-------------------------------------------------")

solve_genetics_cross()