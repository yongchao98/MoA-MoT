import itertools
from collections import Counter
import math

def get_phenotype(genotype):
    """
    Determines the eye color phenotype based on the genotype.
    A fly has wild-type eyes only if it is homozygous recessive (su-v/su-v)
    for the suppressor gene. Otherwise, it is vermilion.
    """
    # Alleles for the suppressor gene are at indices 2 and 3
    autosomal_alleles = (genotype[2], genotype[3])
    if autosomal_alleles.count('su-v') == 2:
        return "wild-type"
    else:
        return "vermilion"

def solve_genetics_cross():
    """
    Solves the Drosophila genetics problem to find the F2 phenotypic ratio.
    """
    # Step 1: Define F1 genotypes based on the problem description.
    # The parental cross is X(v)X(v); su-v/su-v  x  X(v)Y; su-v+/su-v+.
    # This results in an F1 generation where all flies are heterozygous
    # for the su-v gene and have the vermilion genotype.
    f1_female_genotype = ('Xv', 'Xv', 'su-v+', 'su-v')
    f1_male_genotype = ('Xv', 'Y', 'su-v+', 'su-v')

    print("This script calculates the F2 phenotypic ratio from the following F1 intercross:")
    print(f"F1 Female Genotype: {f1_female_genotype[0]}{f1_female_genotype[1]}; {f1_female_genotype[2]}/{f1_female_genotype[3]}")
    print(f"F1 Male Genotype:   {f1_male_genotype[0]}{f1_male_genotype[1]}; {f1_male_genotype[2]}/{f1_male_genotype[3]}\n")

    # Step 2: Generate all possible unique gametes from F1 parents.
    # F1 Female (XvXv; su-v+/su-v) produces two types of gametes.
    female_gametes = [('Xv', allele) for allele in sorted(list(set(f1_female_genotype[2:])))]

    # F1 Male (XvY; su-v+/su-v) produces four types of gametes.
    male_gametes = list(itertools.product(
        sorted(list(set(f1_male_genotype[:2]))),
        sorted(list(set(f1_male_genotype[2:])))
    ))

    print("Possible F1 Female Gametes:", female_gametes)
    print("Possible F1 Male Gametes:  ", male_gametes)
    print("-" * 50)

    # Step 3 & 4: Generate all F2 genotypes, determine phenotypes, and count them.
    phenotype_counts = Counter()

    # The Punnett square is represented by the product of the gamete lists.
    # Each combination represents an equally likely outcome.
    for female_gamete in female_gametes:
        for male_gamete in male_gametes:
            # Combine gametes to form the zygote's genotype.
            sex_alleles = sorted([female_gamete[0], male_gamete[0]])
            if 'Y' in sex_alleles:
                sex_alleles = ['Xv', 'Y']  # Standard representation
            
            autosomal_alleles = sorted([female_gamete[1], male_gamete[1]])
            
            zygote_genotype = tuple(sex_alleles + autosomal_alleles)
            
            # Determine the resulting phenotype and add it to the counter.
            phenotype = get_phenotype(zygote_genotype)
            phenotype_counts[phenotype] += 1
            
    # Step 5: Tally the results.
    vermilion_count = phenotype_counts['vermilion']
    wild_type_count = phenotype_counts['wild-type']
    total_offspring = vermilion_count + wild_type_count

    print("F2 Generation Phenotype Counts:")
    print(f"Total vermilion offspring: {vermilion_count}")
    print(f"Total wild-type offspring: {wild_type_count}")
    print(f"Total possible outcomes: {total_offspring}")
    print("-" * 50)

    # Step 6: Express the result as a simplified ratio and print the final equation.
    common_divisor = math.gcd(vermilion_count, wild_type_count)
    vermilion_ratio_part = vermilion_count // common_divisor
    wild_type_ratio_part = wild_type_count // common_divisor
    total_ratio_parts = total_offspring // common_divisor

    print("Final F2 Phenotypic Ratio Calculation:")
    print(f"The ratio is {vermilion_count} vermilion : {wild_type_count} wild-type, which simplifies to {vermilion_ratio_part}:{wild_type_ratio_part}.")
    
    print("\nFinal Equation:")
    print(f"{vermilion_ratio_part}/{total_ratio_parts} vermilion : {wild_type_ratio_part}/{total_ratio_parts} wild-type")

solve_genetics_cross()