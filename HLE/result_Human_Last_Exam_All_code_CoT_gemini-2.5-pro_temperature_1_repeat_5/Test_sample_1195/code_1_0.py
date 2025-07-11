import math

def solve_genetics_problem():
    """
    Calculates the F2 phenotypic ratio for a cross involving X-linked vermilion
    and autosomal suppressor of vermilion in Drosophila.
    """

    # Step 1: Define the gametes produced by the F1 generation and their probabilities.
    # F1 Female (XvXv; +/su-v) produces two types of gametes.
    f1_female_gametes = {
        ('Xv', '+'): 0.5,      # Gamete with Xv and wild-type suppressor allele
        ('Xv', 'su-v'): 0.5   # Gamete with Xv and recessive suppressor allele
    }

    # F1 Male (XvY; +/su-v) produces four types of gametes.
    f1_male_gametes = {
        ('Xv', '+'): 0.25,     # Gamete with Xv and wild-type suppressor allele
        ('Xv', 'su-v'): 0.25,  # Gamete with Xv and recessive suppressor allele
        ('Y', '+'): 0.25,      # Gamete with Y and wild-type suppressor allele
        ('Y', 'su-v'): 0.25   # Gamete with Y and recessive suppressor allele
    }

    # Step 2: Initialize counters for the F2 phenotypes.
    phenotype_proportions = {
        'wild_type': 0.0,
        'vermilion': 0.0
    }

    # Step 3: Simulate the F1 cross by combining all possible gametes.
    for f_gamete, f_prob in f1_female_gametes.items():
        for m_gamete, m_prob in f1_male_gametes.items():
            # The probability of this specific offspring genotype
            offspring_prob = f_prob * m_prob
            
            # Autosomal alleles from each parent gamete
            autosome_1 = f_gamete[1]
            autosome_2 = m_gamete[1]

            # Determine the phenotype based on the suppressor gene.
            # Suppression occurs (restoring wild-type eyes) only if the genotype is su-v/su-v.
            # All F2 offspring are genetically vermilion (XvXv or XvY), so their
            # final eye color depends entirely on the suppressor gene's action.
            if autosome_1 == 'su-v' and autosome_2 == 'su-v':
                phenotype_proportions['wild_type'] += offspring_prob
            else:
                phenotype_proportions['vermilion'] += offspring_prob

    # Step 4: Convert the final proportions to a simple fraction.
    # The least common multiple of the denominators of probabilities (2 and 4) is 4.
    denominator = 4
    vermilion_numerator = int(phenotype_proportions['vermilion'] * denominator)
    wild_type_numerator = int(phenotype_proportions['wild_type'] * denominator)

    # Step 5: Print the final calculated phenotypic ratio.
    print("The expected phenotypic ratio in the F2 generation is:")
    print(f"{vermilion_numerator}/{denominator} vermilion : {wild_type_numerator}/{denominator} wild-type")

solve_genetics_problem()