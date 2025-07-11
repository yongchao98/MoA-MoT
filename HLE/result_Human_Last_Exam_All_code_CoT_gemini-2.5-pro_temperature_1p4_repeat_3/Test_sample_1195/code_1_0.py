def calculate_drosophila_cross():
    """
    Calculates the F2 phenotypic ratio for a cross involving an X-linked gene
    and an autosomal suppressor gene in Drosophila.
    """
    print("### Solving the Drosophila Genetics Problem Step-by-Step ###\n")

    # Step 1: Parental (P) Generation Genotypes
    print("Step 1: Define Parental (P) Generation Genotypes")
    p_female = "X(v)X(v); su-v/su-v"
    p_male = "X(v)Y; su-v+/su-v+"
    print(f"  - P Female: {p_female} (vermilion female homozygous for su-v)")
    print(f"  - P Male:   {p_male} (vermilion male with wild-type su-v alleles)\n")

    # Step 2: Determine F1 Generation Genotypes
    print("Step 2: Determine F1 Generation Genotypes")
    print("  - All F1 flies will inherit X(v) and su-v from the mother.")
    print("  - All F1 flies will inherit su-v+ from the father.")
    print("  - F1 Females inherit X(v) from the father, F1 Males inherit Y.")
    f1_female = "X(v)X(v); su-v+/su-v"
    f1_male = "X(v)Y; su-v+/su-v"
    print(f"  - F1 Female Genotype: {f1_female}")
    print(f"  - F1 Male Genotype:   {f1_male}")
    print("  - Phenotype: Because all F1 flies are heterozygous (su-v+/su-v), the suppressor is inactive. Therefore, all F1 flies are vermilion.\n")

    # Step 3: Analyze the F2 Generation Cross (F1 x F1)
    print("Step 3: Analyze the F2 Generation Cross")
    print("  - The F2 generation comes from crossing the F1 flies.")
    print("  - We can analyze the X-linked and autosomal genes separately.\n")

    print("  Part A: X-linked cross [X(v)X(v) x X(v)Y]")
    print("  - All F2 offspring will have at least one X(v) allele and no X(+) allele.")
    print("  - Therefore, without suppression, all F2 flies would be vermilion.")
    print("  - This means the final phenotype depends ONLY on the suppressor gene.\n")

    print("  Part B: Autosomal cross [su-v+/su-v x su-v+/su-v]")
    # Standard monohybrid cross ratios
    # Genotypes: 1 su-v+/su-v+, 2 su-v+/su-v, 1 su-v/su-v
    total_outcomes = 4
    num_with_suppression = 1  # su-v/su-v genotype
    num_without_suppression = 3 # su-v+/su-v+ and su-v+/su-v genotypes

    print(f"  - The genotypic ratio is {num_with_suppression+num_without_suppression-1-1}/{total_outcomes} su-v+/su-v+ : {num_with_suppression+num_without_suppression-2}/{total_outcomes} su-v+/su-v : {num_with_suppression}/{total_outcomes} su-v/su-v.\n")

    # Step 4: Calculate the Final F2 Phenotypic Ratio
    print("Step 4: Calculate the Final F2 Phenotypic Ratio")
    print("  - Wild-type eyes appear ONLY when the vermilion color is suppressed (genotype su-v/su-v).")
    print("  - Vermilion eyes appear when there is no suppression (genotypes su-v+/su-v+ or su-v+/su-v).")
    
    wild_type_numerator = num_with_suppression
    vermilion_numerator = num_without_suppression
    denominator = total_outcomes
    
    print("\n--- Final Answer ---")
    print("The expected phenotypic ratio in the F2 generation is:")
    print(f"{vermilion_numerator}/{denominator} vermilion : {wild_type_numerator}/{denominator} wild-type")

# Execute the function to print the result
calculate_drosophila_cross()