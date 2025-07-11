from fractions import Fraction

def get_phenotype(x_geno, auto_geno):
    """Determines the phenotype based on the genotype."""
    # The problem setup ensures all F2 individuals have the 'v' allele.
    # The phenotype is determined solely by the suppressor gene 'su-v'.
    # 'su-v/su-v' is the only genotype that suppresses vermilion, resulting in wild-type eyes.
    if 'su-v' in auto_geno and auto_geno.count('su-v') == 2:
        return "wild-type"
    else:
        return "vermilion"

# F1 Generation Genotypes
f1_female_geno = {'x': ('v', 'v'), 'auto': ('su-v+', 'su-v')}
f1_male_geno = {'x': ('v', 'Y'), 'auto': ('su-v+', 'su-v')}

# Gametes from F1 female (XᵛXᵛ su-v⁺/su-v)
# All X are Xᵛ. Autosomes are 1/2 su-v⁺ and 1/2 su-v
f1_female_gametes = {
    ('v', 'su-v+'): Fraction(1, 2),
    ('v', 'su-v'): Fraction(1, 2)
}

# Gametes from F1 male (XᵛY su-v⁺/su-v)
# 1/2 X are Xᵛ, 1/2 are Y. 1/2 autosomes are su-v⁺, 1/2 are su-v
# Independent assortment gives 4 combinations, each with 1/4 probability.
f1_male_gametes = {
    ('v', 'su-v+'): Fraction(1, 4),
    ('v', 'su-v'): Fraction(1, 4),
    ('Y', 'su-v+'): Fraction(1, 4),
    ('Y', 'su-v'): Fraction(1, 4)
}

# Calculate F2 generation
phenotype_ratios = {"wild-type": Fraction(0), "vermilion": Fraction(0)}

print("Calculating F2 phenotypic ratio from F1 cross...")
print(f"F1 Female (XᵛXᵛ su-v⁺/su-v) x F1 Male (XᵛY su-v⁺/su-v)\n")

for fem_gamete, fem_prob in f1_female_gametes.items():
    for male_gamete, male_prob in f1_male_gametes.items():
        
        # Combine gametes to form F2 offspring genotype
        f2_prob = fem_prob * male_prob
        
        # We only need the autosomal genotype for phenotype determination in this specific problem
        f2_auto_geno = sorted((fem_gamete[1], male_gamete[1]))
        
        # Determine phenotype
        phenotype = get_phenotype(None, f2_auto_geno)
        
        # Add probability to the corresponding phenotype
        phenotype_ratios[phenotype] += f2_prob

# Print the final results
print("Final F2 Phenotypic Ratio:")
total = sum(phenotype_ratios.values())
wt_ratio = phenotype_ratios["wild-type"]
ver_ratio = phenotype_ratios["vermilion"]

print(f"Wild-type: {wt_ratio.numerator}/{wt_ratio.denominator}")
print(f"Vermilion: {ver_ratio.numerator}/{ver_ratio.denominator}")

print(f"\nThe final expected phenotypic ratio is {ver_ratio.numerator}/{total.denominator} vermilion : {wt_ratio.numerator}/{total.denominator} wild-type.")
