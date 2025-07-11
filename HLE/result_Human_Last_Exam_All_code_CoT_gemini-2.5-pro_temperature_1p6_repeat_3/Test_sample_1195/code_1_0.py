def solve_genetics_cross():
    """
    Calculates the F2 phenotypic ratio for the given Drosophila cross.

    The problem involves an X-linked gene (vermilion, v) and an autosomal
    suppressor gene (su-v).

    P Cross: vermilion female (X^vX^v; su-v/su-v) x vermilion male (X^vY; +/+).
    F1 Generation: All offspring are heterozygous for the suppressor (+/su-v) and
                   have the vermilion allele(s), resulting in a vermilion phenotype.
    F1 Cross: X^vX^v; +/su-v x X^vY; +/su-v.

    The F2 phenotype depends only on the segregation of the suppressor gene, as all F2
    flies will inherit the vermilion eye color gene from the F1 parents.
    """

    # We analyze the cross of the heterozygous suppressor gene: +/su-v x +/su-v
    # This gives a 1:2:1 genotypic ratio.
    # Genotypes: 1/4 (+/+), 2/4 (+/su-v), 1/4 (su-v/su-v)

    # Phenotypes:
    # - If 'su-v' is homozygous (su-v/su-v), the suppressor is active,
    #   and the eye color is wild-type.
    # - If at least one dominant '+' allele is present (+/+ or +/su-v),
    #   the suppressor is inactive, and the eye color is vermilion.

    # Fraction of F2 with wild-type eyes (genotype su-v/su-v)
    wild_type_numerator = 1
    denominator = 4

    # Fraction of F2 with vermilion eyes (genotypes +/+ and +/su-v)
    # This is (1/4) + (2/4) = 3/4
    vermilion_numerator = 3

    print("Expected phenotypic ratio in the F2 generation:")
    print(f"Wild-type: {wild_type_numerator} / {denominator}")
    print(f"Vermilion: {vermilion_numerator} / {denominator}")
    print(f"\nThe ratio is {vermilion_numerator}/4 vermilion : {wild_type_numerator}/4 wild-type.")


solve_genetics_cross()
<<<B>>>