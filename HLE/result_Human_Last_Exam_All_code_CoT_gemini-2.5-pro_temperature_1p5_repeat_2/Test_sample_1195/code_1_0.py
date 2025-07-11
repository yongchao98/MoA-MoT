def calculate_f2_phenotypes():
    """
    Calculates the F2 phenotypic ratio for the given Drosophila cross.
    This script breaks down the genetic cross step-by-step.
    """
    print("Step 1: Define Parental (P) Genotypes")
    p_female_genotype = "X(v)X(v); su-v/su-v"
    p_male_genotype = "X(v)Y; su-v+/su-v+"
    print(f"Parental Female: {p_female_genotype}")
    print(f"Parental Male:   {p_male_genotype}\n")

    print("Step 2: Determine F1 Generation Genotypes")
    f1_female_genotype = "X(v)X(v); su-v+/su-v"
    f1_male_genotype = "X(v)Y; su-v+/su-v"
    print(f"All F1 Females are: {f1_female_genotype}")
    print(f"All F1 Males are:   {f1_male_genotype}")
    print("F1 Phenotype: Since the suppressor allele (su-v) is recessive, all F1 flies have vermilion eyes.\n")

    print("Step 3: Analyze the F2 Generation (F1 x F1 Cross)")
    print("The F2 phenotypes depend on the inheritance of the autosomal suppressor gene.")
    print("The cross for this gene is: su-v+/su-v x su-v+/su-v\n")

    # The result of a monohybrid cross su-v+/su-v x su-v+/su-v is:
    # 1/4 su-v+/su-v+ (non-suppressed)
    # 2/4 su-v+/su-v  (non-suppressed)
    # 1/4 su-v/su-v   (suppressed)
    
    # Calculate proportions
    prop_vermilion = 1/4 + 2/4  # su-v+/_ genotypes
    prop_wild_type = 1/4         # su-v/su-v genotype

    # The ratio is 3:1, let's represent the numbers
    vermilion_parts = 3
    wild_type_parts = 1

    print("Step 4: Calculate the final F2 phenotypic ratio.")
    print("All F2 offspring have the vermilion genotype (XvXv or XvY).")
    print("The eye color is determined by the suppressor gene:")
    print(f"- Offspring with su-v/su-v (1/4) are suppressed, resulting in wild-type eyes.")
    print(f"- Offspring with su-v+/_ (3/4) are not suppressed, resulting in vermilion eyes.\n")
    
    print("--- Final F2 Phenotypic Ratio ---")
    print(f"The ratio is {vermilion_parts} vermilion to {wild_type_parts} wild-type.")
    print(f"In fractions, this is {prop_vermilion} vermilion : {prop_wild_type} wild-type.")
    print("\nFinal Equation Numbers:")
    print(f"Vermilion Parts: {vermilion_parts}")
    print(f"Wild-Type Parts: {wild_type_parts}")


# Run the calculation
calculate_f2_phenotypes()