import math
from collections import Counter

def get_phenotype(x_linked_geno, autosomal_geno):
    """
    Determines the eye color phenotype based on the genotype.
    - 'ss' genotype suppresses the vermilion trait, resulting in wild-type eyes.
    - Otherwise, the X-linked gene determines the phenotype.
    """
    if 's' in autosomal_geno and autosomal_geno.count('s') == 2:
        return "wild-type"
    elif 'v' in x_linked_geno:
        return "vermilion"
    else:
        return "wild-type"

# Step 1: Define Parental (P) generation and determine F1 generation
print("--- Step 1: P-generation Cross and F1 Results ---")
p_female_geno = ('XvXv', 'ss')
p_male_geno = ('XvY', 'SS')
print(f"Parental Cross: Female ({p_female_geno[0]}; {p_female_geno[1]}) x Male ({p_male_geno[0]}; {p_male_geno[1]})")

# F1 Generation
f1_female_geno = ('XvXv', 'Ss')
f1_male_geno = ('XvY', 'Ss')
f1_female_pheno = get_phenotype(f1_female_geno[0], f1_female_geno[1])
f1_male_pheno = get_phenotype(f1_male_geno[0], f1_male_geno[1])

print(f"All F1 Females are: Genotype {f1_female_geno[0]}; {f1_female_geno[1]}, Phenotype {f1_female_pheno}")
print(f"All F1 Males are: Genotype {f1_male_geno[0]}; {f1_male_geno[1]}, Phenotype {f1_male_pheno}")
print("\n")


# Step 2: F1 intercross and determine gametes
print("--- Step 2: F1-generation Intercross and Gametes ---")
print(f"F1 Intercross: Female ({f1_female_geno[0]}; {f1_female_geno[1]}) x Male ({f1_male_geno[0]}; {f1_male_geno[1]})")

# Gametes from F1 parents
f1_female_gametes = [('Xv', 'S'), ('Xv', 's')]
f1_male_gametes = [('Xv', 'S'), ('Xv', 's'), ('Y', 'S'), ('Y', 's')]
print(f"F1 Female produces gametes: {f1_female_gametes}")
print(f"F1 Male produces gametes: {f1_male_gametes}")
print("\n")


# Step 3: Determine F2 generation genotypes and phenotypes
print("--- Step 3: F2-generation Genotypes and Phenotypes ---")
f2_phenotypes = []
for f_gamete in f1_female_gametes:
    for m_gamete in f1_male_gametes:
        # Combine gametes to form F2 genotype
        f2_x_geno = f_gamete[0] + m_gamete[0]
        # Sort autosomal alleles for consistency (e.g., 'sS' becomes 'Ss')
        f2_auto_geno = "".join(sorted(f_gamete[1] + m_gamete[1]))
        
        # Determine phenotype
        phenotype = get_phenotype(f2_x_geno, f2_auto_geno)
        f2_phenotypes.append(phenotype)
        print(f"F2 Offspring: Genotype ({f2_x_geno}; {f2_auto_geno}), Phenotype: {phenotype}")

print("\n")

# Step 4: Calculate the final F2 phenotypic ratio
print("--- Step 4: Final F2 Phenotypic Ratio ---")
phenotype_counts = Counter(f2_phenotypes)
total_offspring = len(f2_phenotypes)
gcd = math.gcd(phenotype_counts['wild-type'], phenotype_counts['vermilion'])

wild_type_numerator = phenotype_counts['wild-type'] // gcd
vermilion_numerator = phenotype_counts['vermilion'] // gcd
denominator = total_offspring // gcd

print(f"Counts: {phenotype_counts['wild-type']} wild-type, {phenotype_counts['vermilion']} vermilion")
print(f"Total Offspring: {total_offspring}")
print("\nThe expected phenotypic ratio in the F2 generation is:")
print(f"{vermilion_numerator}/{denominator} vermilion : {wild_type_numerator}/{denominator} wild-type")
<<<B>>>