def solve_drosophila_cross():
    """
    Calculates and explains the F2 phenotypic ratio for the given Drosophila cross.
    """
    print("This script calculates the phenotypic ratio for a dihybrid cross in Drosophila.")
    print("The cross involves an X-linked gene (vermilion) and an autosomal suppressor gene (su-v).\n")
    print("Interaction Rule: The 'su-v/su-v' genotype restores wild-type eye color, even in flies with the vermilion genotype.")

    # Step 1: Define Parental Generation
    print("--- Step 1: Parental (P) Generation ---")
    p_female_genotype = "X(v)X(v) ; su-v/su-v"
    p_male_genotype = "X(v)Y ; +/+"
    print(f"Parental Female Genotype: {p_female_genotype}")
    print(f"Parental Male Genotype:   {p_male_genotype}")
    print("(Note: Based on the rules, a female with this genotype should be wild-type. We proceed assuming the provided genotypes for the cross are correct.)")

    # Step 2: Determine F1 Generation
    print("\n--- Step 2: F1 Generation (from P cross) ---")
    f1_female_genotype = "X(v)X(v) ; +/su-v"
    f1_male_genotype = "X(v)Y ; +/su-v"
    print(f"All F1 Females are: {f1_female_genotype}")
    print(f"All F1 Males are:   {f1_male_genotype}")
    print("Phenotype: All F1 flies are vermilion because none have the 'su-v/su-v' genotype required for suppression.")

    # Step 3: F1 Cross and F2 Analysis
    print("\n--- Step 3: F2 Generation (from F1 x F1 cross) ---")
    print("The F1 generation is intercrossed. We analyze the inheritance of each gene separately.")
    
    print("\n- Analysis of X-linked gene (X(v)X(v) x X(v)Y):")
    print("  All F2 offspring have the genetic basis for vermilion eyes (X(v)X(v) or X(v)Y).")
    
    print("\n- Analysis of autosomal suppressor gene (+/su-v x +/su-v):")
    # A standard monohybrid cross results in a 1:2:1 genotypic ratio.
    genotypes_su_v = {
        "+/+": 1,
        "+/su-v": 2,
        "su-v/su-v": 1
    }
    total_parts = sum(genotypes_su_v.values())
    print(f"  The F2 genotypic ratio for the suppressor gene is: {genotypes_su_v['+/+']} (+/+) : {genotypes_su_v['+/su-v']} (+/su-v) : {genotypes_su_v['su-v/su-v']} (su-v/su-v).")

    # Step 4: Calculate Final F2 Phenotypic Ratio
    print("\n--- Step 4: Final F2 Phenotypic Ratio ---")
    print("The final eye color depends on whether the vermilion trait is suppressed:")
    
    # Calculate phenotypes based on the epistatic interaction
    # Vermilion phenotype: occurs in +/+ and +/su-v genotypes
    vermilion_parts = genotypes_su_v["+/+"] + genotypes_su_v["+/su-v"]
    
    # Wild-type phenotype: occurs in su-v/su-v genotype (suppression)
    wild_type_parts = genotypes_su_v["su-v/su-v"]
    
    print(f"- Individuals with genotype 'su-v/su-v' have their vermilion color suppressed, resulting in wild-type eyes.")
    print(f"- Individuals with genotypes '+/+' or '+/su-v' are not suppressed, resulting in vermilion eyes.")

    print("\nFinal Equation:")
    print(f"The total ratio is composed of {total_parts} parts.")
    print(f"Number of vermilion parts = {vermilion_parts}")
    print(f"Number of wild-type parts = {wild_type_parts}")
    
    print(f"\nThe expected phenotypic ratio in the F2 generation is:")
    print(f"{vermilion_parts} vermilion : {wild_type_parts} wild-type")
    print(f"This simplifies to a fractional ratio of {vermilion_parts}/{total_parts} vermilion to {wild_type_parts}/{total_parts} wild-type.")

solve_drosophila_cross()
<<<B>>>