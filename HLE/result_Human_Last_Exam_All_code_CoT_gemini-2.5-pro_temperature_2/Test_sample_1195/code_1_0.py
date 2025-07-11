import math

def solve_drosophila_cross():
    """
    Calculates the F2 phenotypic ratio for a cross involving an X-linked gene
    and an autosomal suppressor gene.
    """

    # In the F1 x F1 cross, the determining factor for the phenotype is the
    # monohybrid cross of the suppressor gene: su-v+/su-v  x  su-v+/su-v
    # The outcome of this cross follows a 1:2:1 genotypic ratio.
    # Total parts in the ratio = 1 + 2 + 1 = 4.
    total_outcomes = 4

    # 1 part of the offspring will have the genotype su-v+/su-v+
    # 2 parts will have the genotype su-v+/su-v
    # 1 part will have the genotype su-v/su-v

    # Phenotype determination:
    # - If the genotype is su-v/su-v, it suppresses the vermilion allele,
    #   and the fly has wild-type eyes.
    wild_type_count = 1

    # - If the genotype has at least one su-v+ allele (su-v+/su-v+ or su-v+/su-v),
    #   suppression does not occur, and the fly has vermilion eyes.
    vermilion_count = 1 + 2

    # Simplify the fractions to their lowest common denominator
    common_divisor = math.gcd(wild_type_count, vermilion_count)
    if common_divisor > 1: # Though in this case it's 1
        wild_type_num = wild_type_count // common_divisor
        vermilion_num = vermilion_count // common_divisor
    else:
        wild_type_num = wild_type_count
        vermilion_num = vermilion_count

    print("The cross for the suppressor gene is su-v+/su-v x su-v+/su-v.")
    print(f"The resulting genotypic ratio is 1 su-v+/su-v+ : 2 su-v+/su-v : 1 su-v/su-v.")
    print("Phenotype Mapping:")
    print("- Genotypes su-v+/su-v+ and su-v+/su-v do not suppress vermilion.")
    print("- Genotype su-v/su-v suppresses vermilion, resulting in wild-type eyes.")
    print("\nExpected F2 phenotypic ratio:")
    # The final equation showing the number of parts for each phenotype
    print(f"{wild_type_count}/{total_outcomes} wild-type : {vermilion_count}/{total_outcomes} vermilion")
    # Also expressed as a simplified ratio
    print(f"Which simplifies to {wild_type_num} wild-type : {vermilion_num} vermilion.")

solve_drosophila_cross()