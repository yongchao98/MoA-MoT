import math

def solve_genetics_cross():
    """
    Calculates the F2 phenotypic ratio for the given Drosophila cross.
    """
    # Step 1: Define the gametes from the F1 generation.
    # F1 Female (XvXv; su-v+/su-v) produces two gamete types:
    f1_female_gametes = [('Xv', 'su-v+'), ('Xv', 'su-v')]

    # F1 Male (XvY; su-v+/su-v) produces four gamete types:
    f1_male_gametes = [('Xv', 'su-v+'), ('Xv', 'su-v'), ('Y', 'su-v+'), ('Y', 'su-v')]

    # Step 2: Initialize counters for the F2 phenotypes.
    phenotypes = {
        'wild_type': 0,
        'vermilion': 0
    }

    # Step 3: Simulate the F1 cross (Punnett square) to get F2 outcomes.
    for female_gamete in f1_female_gametes:
        for male_gamete in f1_male_gametes:
            # All offspring have the vermilion genotype (XvXv or XvY).
            # The phenotype is determined by the su-v gene.
            
            # Combine the su-v alleles from the gametes.
            su_v_allele1 = female_gamete[1]
            su_v_allele2 = male_gamete[1]
            
            # Check for homozygous recessive su-v, which suppresses vermilion.
            if su_v_allele1 == 'su-v' and su_v_allele2 == 'su-v':
                # Suppression occurs, resulting in wild-type eyes.
                phenotypes['wild_type'] += 1
            else:
                # No suppression, resulting in vermilion eyes.
                phenotypes['vermilion'] += 1

    # Step 4: Calculate the simplified ratio.
    total_outcomes = len(f1_female_gametes) * len(f1_male_gametes)
    wild_type_count = phenotypes['wild_type']
    vermilion_count = phenotypes['vermilion']
    
    # Find the greatest common divisor to simplify the ratio.
    common_divisor = math.gcd(wild_type_count, vermilion_count)
    
    wild_type_ratio = wild_type_count // common_divisor
    vermilion_ratio = vermilion_count // common_divisor

    # Step 5: Print the final result.
    print(f"The analysis of the F1 cross results in {total_outcomes} possible F2 genotypes.")
    print(f"Out of these possibilities:")
    print(f"- {wild_type_count} result in wild-type eyes (due to suppression).")
    print(f"- {vermilion_count} result in vermilion eyes (no suppression).")
    print("\nThe final expected phenotypic ratio is:")
    print(f"{vermilion_ratio} vermilion : {wild_type_ratio} wild-type")
    print(f"\nThis corresponds to a ratio of 3/4 vermilion to 1/4 wild-type.")

solve_genetics_cross()