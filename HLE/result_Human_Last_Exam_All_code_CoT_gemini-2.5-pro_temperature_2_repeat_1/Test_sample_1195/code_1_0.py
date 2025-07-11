def solve_genetics_problem():
    """
    Calculates and explains the F2 phenotypic ratio for a Drosophila cross
    involving X-linked vermilion (v) and autosomal suppressor of vermilion (su-v).
    """

    print("Step 1: Define Parental (P) Generation and F1 Offspring")
    print("-" * 60)
    print("Parental Cross:")
    print("  Female Genotype: X(v)X(v) ; su-v/su-v  (Phenotype: Wild-type due to suppression)")
    print("  Male Genotype:   X(v)Y ; su-v+/su-v+ (Phenotype: Vermilion)")
    print()

    # Determine F1 Generation
    f1_female_genotype = "X(v)X(v) ; su-v+/su-v"
    f1_male_genotype = "X(v)Y ; su-v+/su-v"

    print("F1 Generation Offspring:")
    print(f"  All F1 Females: {f1_female_genotype}")
    print(f"  All F1 Males:   {f1_male_genotype}")
    print("  F1 Phenotype: Since all F1 flies are heterozygous for the suppressor gene (su-v+/su-v),")
    print("  the recessive suppressor is not expressed. All F1 flies show the vermilion eye color.")
    print()

    print("Step 2: Analyze the F1 Intercross to determine the F2 Generation")
    print("-" * 60)
    print(f"F1 Intercross: (Female) {f1_female_genotype} x (Male) {f1_male_genotype}")
    print()
    print("In this cross, all offspring will have the genotype for vermilion eyes (X(v)X(v) or X(v)Y).")
    print("The final eye color phenotype depends only on the segregation of the autosomal suppressor gene.")
    print()
    print("The cross for the suppressor gene is a standard monohybrid cross: su-v+/su-v x su-v+/su-v")
    print("This yields the following genotypic ratio in the F2 generation:")

    # Standard Mendelian ratio for the autosomal gene
    f2_auto_genotypes = {
        "su-v+/su-v+": 1,
        "su-v+/su-v": 2,
        "su-v/su-v": 1
    }
    total_parts = sum(f2_auto_genotypes.values())

    print(f"  {f2_auto_genotypes['su-v+/su-v+']}/{total_parts} su-v+/su-v+")
    print(f"  {f2_auto_genotypes['su-v+/su-v']}/{total_parts} su-v+/su-v")
    print(f"  {f2_auto_genotypes['su-v/su-v']}/{total_parts} su-v/su-v")
    print()

    print("Step 3: Calculate the F2 Phenotypic Ratio")
    print("-" * 60)
    # Phenotype determination
    # Genotypes with at least one su-v+ do not suppress -> Vermilion
    # Genotype su-v/su-v suppresses -> Wild-type
    vermilion_count = f2_auto_genotypes["su-v+/su-v+"] + f2_auto_genotypes["su-v+/su-v"]
    wild_type_count = f2_auto_genotypes["su-v/su-v"]

    print(f"{vermilion_count}/{total_parts} of offspring have at least one dominant 'su-v+' allele.")
    print("  - Result: No suppression occurs. Phenotype is Vermilion.")
    print()
    print(f"{wild_type_count}/{total_parts} of offspring are homozygous recessive 'su-v/su-v'.")
    print("  - Result: Suppression occurs. Phenotype is Wild-type.")
    print()

    print("Final F2 Phenotypic Ratio:")
    print("The final equation for the ratio of eye color phenotypes is:")
    print(f"{vermilion_count}/{total_parts} vermilion : {wild_type_count}/{total_parts} wild-type")

# Execute the function to print the solution
solve_genetics_problem()