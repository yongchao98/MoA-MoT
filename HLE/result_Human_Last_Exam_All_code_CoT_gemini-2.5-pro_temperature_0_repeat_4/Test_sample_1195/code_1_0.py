def solve_genetics_problem():
    """
    This script solves the Drosophila eye color genetics problem by calculating the F2 phenotypic ratio.

    The plan is as follows:
    1.  Acknowledge and resolve the contradiction in the problem's description of the parental female.
        The most logical cross is a wild-type female (XᵛXᵛ; su-v/su-v) with a vermilion male (XᵛY; su-v+/su-v+).
    2.  Determine the F1 generation's genotype. All F1s are heterozygous for the suppressor gene (su-v+/su-v)
        and have the vermilion allele, making them phenotypically vermilion.
    3.  Analyze the F2 generation from an F1 intercross. Since all F2s have the vermilion allele (Xᵛ),
        the phenotype depends only on the monohybrid cross of the suppressor gene: su-v+/su-v x su-v+/su-v.
    4.  Calculate the resulting phenotypic ratio based on the principles of recessive epistasis.
    """

    # The F2 generation's phenotype depends on the outcome of the su-v+/su-v x su-v+/su-v cross.
    # The genotypic ratio is 1 (su-v+/su-v+) : 2 (su-v+/su-v) : 1 (su-v/su-v).
    # Total parts in the ratio is 1 + 2 + 1 = 4.
    total_parts = 4

    # Count genotypes that result in the vermilion phenotype (no suppression).
    # These are su-v+/su-v+ (1 part) and su-v+/su-v (2 parts).
    vermilion_parts = 1 + 2

    # Count genotypes that result in the wild-type phenotype (suppression).
    # This is su-v/su-v (1 part).
    wild_type_parts = 1

    # Print the final equation representing the phenotypic ratio.
    print("The expected phenotypic ratio in the F2 generation is:")
    print(f"{vermilion_parts}/{total_parts} vermilion : {wild_type_parts}/{total_parts} wild-type")
    print("\nThe final equation is:")
    print(f"{vermilion_parts} vermilion : {wild_type_parts} wild-type")

solve_genetics_problem()
<<<B>>>