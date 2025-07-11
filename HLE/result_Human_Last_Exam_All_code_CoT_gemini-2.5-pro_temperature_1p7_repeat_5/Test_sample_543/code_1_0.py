# This script determines and prints the IUPAC name of the reaction product.

# The reaction involves 1,3-dibromo-2-iodobenzene with excess phenyl magnesium bromide
# under reflux conditions, followed by aqueous work-up. These forcing conditions
# drive the reaction to substitute all three halogen atoms with phenyl groups.

# The original positions of the halogens are 1, 2, and 3.
# The reaction replaces each halogen with a phenyl group.
parent_molecule = "benzene"
substituent = "phenyl"
count_prefix = "tri"
positions = [1, 2, 3]

# Constructing the IUPAC name
position_str = ",".join(map(str, positions))
final_name = f"{position_str}-{count_prefix}{substituent}{parent_molecule}"

# Print the final result. The numbers 1, 2, and 3 are included in the output.
print("The IUPAC name of the final product is:")
print(final_name)