from collections import defaultdict

def get_phenotype(sex_chromosome, autosomal_genotype):
    """Determines the eye color phenotype based on genotype."""
    # All individuals in the F2 generation will have the vermilion allele (Xv).
    # The phenotype depends solely on the suppressor gene.
    # su-v/su-v genotype suppresses vermilion, resulting in wild-type eyes.
    if 'su-v/su-v' in autosomal_genotype:
        return "wild-type"
    else:
        # su-v+/su-v+ or su-v+/su-v genotypes do not suppress vermilion.
        return "vermilion"

def solve_drosophila_cross():
    """
    Calculates the F2 phenotypic ratio for the given Drosophila cross.
    """
    print("--- F1 Generation Analysis ---")
    # P-gen: Female (XvXv; su-v/su-v) x Male (XvY; su-v+/su-v+)
    # All F1 are heterozygous for the su-v gene.
    f1_female_genotype = ('XvXv', 'su-v+/su-v')
    f1_male_genotype = ('XvY', 'su-v+/su-v')
    print(f"F1 Female Genotype: {f1_female_genotype[0]}; {f1_female_genotype[1]}")
    print(f"F1 Male Genotype:   {f1_male_genotype[0]}; {f1_male_genotype[1]}")
    print("All F1 generation flies have vermilion eyes.\n")


    print("--- F2 Generation Analysis ---")
    # Gametes from F1 female (XvXv; su-v+/su-v)
    # She is homozygous for Xv, so all gametes get Xv.
    # She is heterozygous for su-v, so gametes get su-v+ or su-v.
    f1_female_gametes = {
        ('Xv', 'su-v+'): 0.5,
        ('Xv', 'su-v'): 0.5
    }
    print("F1 Female Gametes and Probabilities:", f1_female_gametes)

    # Gametes from F1 male (XvY; su-v+/su-v)
    f1_male_gametes = {
        ('Xv', 'su-v+'): 0.25,
        ('Xv', 'su-v'): 0.25,
        ('Y', 'su-v+'): 0.25,
        ('Y', 'su-v'): 0.25
    }
    print("F1 Male Gametes and Probabilities:", f1_male_gametes)
    print("\n--- Calculating F2 Phenotypic Ratio ---")

    phenotype_ratios = defaultdict(float)
    total_offspring = 0

    # Simulate the Punnett Square
    for f_gamete, f_prob in f1_female_gametes.items():
        for m_gamete, m_prob in f1_male_gametes.items():
            sex_chromo = f_gamete[0] + m_gamete[0]
            # Sort alleles for consistent representation, e.g., su-v+/su-v
            auto_alleles = sorted([f_gamete[1], m_gamete[1]], reverse=True)
            autosomal_genotype = '/'.join(auto_alleles)

            phenotype = get_phenotype(sex_chromo, autosomal_genotype)
            probability = f_prob * m_prob
            phenotype_ratios[phenotype] += probability

    vermilion_ratio = phenotype_ratios["vermilion"]
    wild_type_ratio = phenotype_ratios["wild-type"]

    # The common denominator for 0.25 is 4.
    denominator = 4
    vermilion_numerator = int(vermilion_ratio * denominator)
    wild_type_numerator = int(wild_type_ratio * denominator)

    print(f"The total proportion of vermilion flies is {vermilion_ratio}")
    print(f"The total proportion of wild-type flies is {wild_type_ratio}\n")

    print("Final F2 Phenotypic Ratio:")
    print(f"{vermilion_numerator}/{denominator} vermilion : {wild_type_numerator}/{denominator} wild-type")


solve_drosophila_cross()
<<<B>>>