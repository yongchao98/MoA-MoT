import collections
from fractions import Fraction

# This script calculates the F2 phenotypic ratio from a Drosophila cross
# involving an X-linked gene (vermilion) and an autosomal suppressor gene (su-v).

# 1. Define the gametes from the F1 parents.
# The F1 generation results from crossing P: (X^vX^v su-v/su-v) x (X^vY su-v+/su-v+).
# This gives F1 genotypes: (Female) X^vX^v su-v+/su-v and (Male) X^vY su-v+/su-v.

# Gametes from F1 Female (X^vX^v su-v+/su-v):
f1_female_gametes = [('Xv', 'su-v+'), ('Xv', 'su-v')]

# Gametes from F1 Male (X^vY su-v+/su-v):
f1_male_gametes = [('Xv', 'su-v+'), ('Xv', 'su-v'), ('Y', 'su-v+'), ('Y', 'su-v')]

# 2. Generate all possible F2 offspring genotypes by combining gametes.
f2_offspring_genotypes = []
for female_gamete in f1_female_gametes:
    for male_gamete in f1_male_gametes:
        f2_offspring_genotypes.append((female_gamete, male_gamete))

# 3. Determine the phenotype for each F2 offspring and count them.
# Phenotype Rule: If the autosomal genotype is 'su-v/su-v', the phenotype is wild-type.
# Otherwise, the vermilion phenotype from the 'Xv' allele is expressed.
phenotype_counts = collections.Counter()
for offspring_gametes in f2_offspring_genotypes:
    autosomal_allele1 = offspring_gametes[0][1]
    autosomal_allele2 = offspring_gametes[1][1]
    
    # Check for homozygous recessive suppressor genotype
    if autosomal_allele1 == 'su-v' and autosomal_allele2 == 'su-v':
        phenotype = 'Wild-type'
    else:
        phenotype = 'Vermilion'
    
    phenotype_counts[phenotype] += 1

# 4. Calculate and print the final phenotypic ratio.
total_offspring = len(f2_offspring_genotypes)
wild_type_count = phenotype_counts['Wild-type']
vermilion_count = phenotype_counts['Vermilion']

print("F2 Generation Phenotypic Ratio Calculation:")
print("-" * 45)

# Output the individual numbers that form the final ratio equation
print(f"Number of possible F2 outcomes: {total_offspring}")
print(f"Number of Vermilion offspring: {vermilion_count}")
print(f"Number of Wild-type offspring: {wild_type_count}")
print("-" * 45)

# Output the final equation
print("Final Ratio Equation:")
wild_type_fraction = Fraction(wild_type_count, total_offspring)
vermilion_fraction = Fraction(vermilion_count, total_offspring)

print(f"Vermilion: {vermilion_count}/{total_offspring} = {vermilion_fraction}")
print(f"Wild-type: {wild_type_count}/{total_offspring} = {wild_type_fraction}")
print("\nThe expected phenotypic ratio is 3/4 vermilion: 1/4 wild-type.")
