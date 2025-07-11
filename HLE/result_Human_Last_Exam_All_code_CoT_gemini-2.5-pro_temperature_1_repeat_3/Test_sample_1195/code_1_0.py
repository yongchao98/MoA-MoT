# This script calculates the F2 phenotypic ratio for the described cross.

# Step 1: Analyze the F1 cross for the autosomal suppressor gene.
# The F1 generation is heterozygous for the suppressor gene (su-v+/su-v).
# The F1 cross is therefore: su-v+/su-v x su-v+/su-v.
# This is a standard monohybrid cross.
# Let's define the parts of the resulting 1:2:1 genotypic ratio.
homozygous_dominant_parts = 1  # su-v+/su-v+
heterozygous_parts = 2         # su-v+/su-v
homozygous_recessive_parts = 1 # su-v/su-v
total_parts = homozygous_dominant_parts + heterozygous_parts + homozygous_recessive_parts

# Step 2: Determine the phenotypes based on the genotypes.
# All F2 flies are genotypically vermilion (XvXv or XvY).
# The phenotype is determined by the suppressor gene.
# - The recessive homozygous genotype 'su-v/su-v' suppresses vermilion, causing wild-type eyes.
# - The 'su-v+' allele is dominant, so 'su-v+/su-v+' and 'su-v+/su-v' genotypes do not suppress, causing vermilion eyes.

wild_type_phenotype_parts = homozygous_recessive_parts
vermilion_phenotype_parts = homozygous_dominant_parts + heterozygous_parts

# Step 3: Print the final phenotypic ratio.
print("F2 Phenotypic Ratio Calculation:")
print("All F2 offspring have a genotype that would result in vermilion eyes (XvXv or XvY).")
print("Their final phenotype depends on the suppressor gene (su-v).")
print(f"The cross for the suppressor gene (su-v+/su-v x su-v+/su-v) yields {total_parts} total parts in the ratio.")
print("-" * 30)
print("Phenotype Breakdown:")
print(f"Number of parts resulting in Wild-type eyes (from 'su-v/su-v' genotype): {wild_type_phenotype_parts}")
print(f"Number of parts resulting in Vermilion eyes (from 'su-v+/su-v+' and 'su-v+/su-v' genotypes): {vermilion_phenotype_parts}")
print("-" * 30)
print("The final expected phenotypic ratio is:")
print(f"{vermilion_phenotype_parts} vermilion : {wild_type_phenotype_parts} wild-type")
print(f"(This is equivalent to 3/4 vermilion : 1/4 wild-type)")