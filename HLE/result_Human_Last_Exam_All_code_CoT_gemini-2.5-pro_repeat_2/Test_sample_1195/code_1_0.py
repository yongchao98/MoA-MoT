# This script calculates the F2 phenotypic ratio for the given Drosophila cross.

# Step 1: Analyze the F1 cross for the autosomal suppressor gene (su-v).
# The F1 generation is a cross between su-v+/su-v and su-v+/su-v.
# We can model the offspring proportions from this monohybrid cross.
# Proportions of genotypes in F2:
# 1 part su-v+/su-v+
# 2 parts su-v+/su-v
# 1 part su-v/su-v
# Total parts = 1 + 2 + 1 = 4

total_parts = 4

# Step 2: Determine the phenotype based on the su-v genotype.
# The X-linked gene is always Xv, so the phenotype is controlled by the suppressor.
# - The su-v/su-v genotype suppresses the vermilion trait, resulting in wild-type eyes.
#   This corresponds to 1 part of the offspring.
wild_type_numerator = 1

# - The su-v+/su-v+ and su-v+/su-v genotypes do not suppress vermilion, resulting in vermilion eyes.
#   This corresponds to 1 + 2 = 3 parts of the offspring.
vermilion_numerator = 3

# Step 3: Print the final phenotypic ratio.
print("The expected F2 phenotypic ratio involves two phenotypes: vermilion and wild-type.")
print("The final equation for the ratio is:")
print(f"{vermilion_numerator} / {total_parts} vermilion : {wild_type_numerator} / {total_parts} wild-type")
