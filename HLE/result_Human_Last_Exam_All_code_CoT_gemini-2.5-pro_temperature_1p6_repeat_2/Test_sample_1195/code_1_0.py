# This script calculates the expected F2 phenotypic ratio based on the provided genetic cross.

# Step 1: Define the results from the F1 monohybrid cross for the suppressor gene (su-v+/su-v x su-v+/su-v).
# The cross for the X-linked gene (XvXv x XvY) results in all offspring being genotypically vermilion.
# Therefore, the final phenotype depends only on the suppressor gene's genotype.
# The total number of parts in the ratio is 4.
denominator = 4

# The genotype su-v/su-v causes suppression, leading to a wild-type phenotype.
# Its probability in a monohybrid cross is 1/4.
wild_type_numerator = 1

# The genotypes su-v+/su-v+ and su-v+/su-v do not cause suppression, leading to a vermilion phenotype.
# Their combined probability is 1/4 + 2/4 = 3/4.
vermilion_numerator = 3

# Step 2: Print the final calculated phenotypic ratio.
# The problem asks for the ratio of vermilion to wild-type.
print("The expected phenotypic ratio in the F2 generation is:")
print(f"{vermilion_numerator}/{denominator} vermilion : {wild_type_numerator}/{denominator} wild-type")
