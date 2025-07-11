# The problem simplifies to a standard monohybrid cross for the suppressor gene (su-v),
# as all F1 individuals are genotypically vermilion (v) and heterozygous
# for the suppressor (su-v+/su-v). The F2 phenotypes depend on the inheritance
# of the su-v gene.

# The F1 cross for the autosomal gene is: su-v+/su-v  x  su-v+/su-v

# In a monohybrid cross, the genotypic ratio of the offspring is 1:2:1.
# Let's define the parts of this ratio.
homozygous_dominant_part = 1  # Represents su-v+/su-v+
heterozygous_part = 2        # Represents su-v+/su-v
homozygous_recessive_part = 1 # Represents su-v/su-v
total_parts = homozygous_dominant_part + heterozygous_part + homozygous_recessive_part

# The phenotype is determined by whether the vermilion color is suppressed.
# Suppression only occurs in the homozygous recessive (su-v/su-v) state.
# All other genotypes will result in the vermilion phenotype.
vermilion_parts = homozygous_dominant_part + heterozygous_part
wild_type_parts = homozygous_recessive_part

print("Calculating the F2 phenotypic ratio based on the cross: su-v+/su-v x su-v+/su-v")
print("----------------------------------------------------------------------")
print("Phenotype Determination:")
print("- Genotypes (su-v+/su-v+) and (su-v+/su-v) are NOT suppressed and result in VERMILION eyes.")
print("- Genotype (su-v/su-v) IS suppressed and results in WILD-TYPE eyes.")
print("\nCalculating the parts for each phenotype from the 1:2:1 genotypic ratio:")
print(f"Parts for vermilion = {homozygous_dominant_part} (from su-v+/su-v+) + {heterozygous_part} (from su-v+/su-v) = {vermilion_parts} parts")
print(f"Parts for wild-type = {homozygous_recessive_part} (from su-v/su-v) = {wild_type_parts} part")
print(f"\nTotal parts in ratio = {total_parts}")
print("\nFinal Equation for Phenotypic Ratio:")
print(f"{vermilion_parts}/{total_parts} vermilion : {wild_type_parts}/{total_parts} wild-type")
print(f"This is equivalent to 3/4 vermilion : 1/4 wild-type.")
