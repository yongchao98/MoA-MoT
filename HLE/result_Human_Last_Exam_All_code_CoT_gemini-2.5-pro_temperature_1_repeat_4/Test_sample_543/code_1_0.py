# This script determines and prints the IUPAC name of the reaction product.

# The reaction involves the substitution of all three halogens (one iodine and
# two bromines) on the 1,3-dibromo-2-iodobenzene molecule with phenyl groups
# from the excess phenyl magnesium bromide.

# The final product has a central benzene ring with phenyl substituents
# at positions 1, 2, and 3. We will now construct the IUPAC name for it.

# Define the components of the IUPAC name
locant_1 = 1
locant_2 = 2
locant_3 = 3
prefix = "tri"
substituent = "phenyl"
parent_ring = "benzene"

# Construct the final IUPAC name by combining the components.
# This method explicitly uses each number as requested.
final_iupac_name = f"{locant_1},{locant_2},{locant_3}-{prefix}{substituent}{parent_ring}"

# Print the final name
print("The IUPAC name of the product is:")
print(final_iupac_name)