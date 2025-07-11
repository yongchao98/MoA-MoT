def solve_genetics_problem():
    """
    This script solves the Drosophila gene interaction problem step-by-step.
    """

    print("Step 1: Define the Parental (P) generation cross.")
    p_female_genotype = "XvXv; su-v/su-v"
    p_male_genotype = "XvY; +/+"
    print(f"Parental Female Genotype: {p_female_genotype}")
    print(f"Parental Male Genotype:   {p_male_genotype}")
    print("-" * 30)

    print("Step 2: Determine the F1 generation.")
    f1_female_genotype = "XvXv; +/su-v"
    f1_male_genotype = "XvY; +/su-v"
    print("All F1 offspring are heterozygous for the suppressor gene (+/su-v).")
    print("Since the su-v allele is recessive, it is not expressed.")
    print("All F1 flies have the vermilion phenotype.")
    print(f"F1 Cross: {f1_female_genotype} (Female) x {f1_male_genotype} (Male)")
    print("-" * 30)

    print("Step 3: Analyze the F2 generation based on the suppressor gene.")
    print("The eye color phenotype in the F2 generation depends solely on the su-v gene,")
    print("as all flies possess the vermilion (Xv) allele.")
    print("We analyze the results of the cross: +/su-v  x  +/su-v")
    print("\nThe resulting F2 genotypic ratio for the su-v gene is 1 (+/+) : 2 (+/su-v) : 1 (su-v/su-v).")
    print("-" * 30)

    print("Step 4: Calculate the final phenotypic ratio.")
    # Genotypes that lead to the Vermilion phenotype
    vermilion_ratio_numerator = 1  # for +/+
    vermilion_ratio_numerator += 2 # for +/su-v
    total_parts = 4

    # Genotypes that lead to the Wild-type phenotype
    wild_type_ratio_numerator = 1 # for su-v/su-v

    print("Phenotype Determination:")
    print(f"- Genotypes +/+ and +/su-v do not suppress. This accounts for {vermilion_ratio_numerator}/{total_parts} of the offspring.")
    print("  These flies will be VERMILION.")
    print(f"- Genotype su-v/su-v suppresses the vermilion trait. This accounts for {wild_type_ratio_numerator}/{total_parts} of the offspring.")
    print("  These flies will be WILD-TYPE.")
    print("-" * 30)
    
    print("Final F2 Phenotypic Ratio Equation:")
    # Using variables to construct the final output string
    vermilion_fraction = f"{vermilion_ratio_numerator}/{total_parts}"
    wild_type_fraction = f"{wild_type_ratio_numerator}/{total_parts}"
    
    print(f"The expected ratio is {vermilion_fraction} vermilion : {wild_type_fraction} wild-type.")
    print("\nFinal Equation:")
    print(f"Total Offspring = {vermilion_fraction} (vermilion) + {wild_type_fraction} (wild-type)")


solve_genetics_problem()
<<<B>>>