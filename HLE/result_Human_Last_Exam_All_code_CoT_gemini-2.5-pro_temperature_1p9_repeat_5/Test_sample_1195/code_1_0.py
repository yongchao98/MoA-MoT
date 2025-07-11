def solve_genetics_problem():
    """
    Calculates and explains the F2 phenotypic ratio for the given Drosophila cross.
    """
    # Define phenotypes
    phenotype_vermilion = "vermilion"
    phenotype_wild_type = "wild-type"

    # The F2 phenotype depends on the segregation of the autosomal su-v gene,
    # as all F2 flies inherit the X-linked vermilion (v) gene.
    # The F1 cross for the autosomal gene is: +/su-v  x  +/su-v

    # This is a standard monohybrid cross with genotypic ratio 1:2:1
    # We can represent the total parts of the ratio as 4.
    
    # Genotypes with at least one dominant '+' allele do not suppress vermilion.
    # These are +/+ (1 part) and +/su-v (2 parts).
    vermilion_parts = 1 + 2

    # Genotypes with homozygous recessive 'su-v/su-v' (1 part) suppress vermilion,
    # resulting in a wild-type phenotype.
    wild_type_parts = 1

    total_parts = vermilion_parts + wild_type_parts

    print("Analysis of the F2 generation based on the F1 cross (+/su-v x +/su-v):")
    print("-" * 60)
    print(f"The F2 phenotypic ratio is determined out of {total_parts} total parts.")
    
    # Final Equation Output
    print(f"\nEquation for Vermilion flies: (1 part for +/+) + (2 parts for +/su-v) = {vermilion_parts} parts")
    print(f"Equation for Wild-Type flies: (1 part for su-v/su-v) = {wild_type_parts} part")

    print("\nResulting F2 Phenotypic Ratio:")
    print(f"{vermilion_parts}/{total_parts} {phenotype_vermilion} : {wild_type_parts}/{total_parts} {phenotype_wild_type}")
    print(f"This simplifies to a ratio of 3/4 {phenotype_vermilion} to 1/4 {phenotype_wild_type}.")

solve_genetics_problem()
<<<B>>>