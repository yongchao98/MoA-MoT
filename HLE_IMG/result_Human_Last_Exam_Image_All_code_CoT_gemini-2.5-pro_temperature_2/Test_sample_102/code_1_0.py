import sys

# This script is designed to output the IUPAC name of the reaction product.
# The reaction is identified as an Oxy-Cope rearrangement.
# The final product is a 10-membered ring ketone.

# Step 1: Define the components of the IUPAC name.
product_name_parts = {
    "stereoisomer": "(E)",
    "substituent_1": "2-methoxy",
    "substituent_2": "10-methyl",
    "parent_ring": "cyclodec",
    "unsaturation": "2-en",
    "principal_group": "1-one",
}

# Step 2: Assemble the full name string.
full_name = f"{product_name_parts['stereoisomer']}-{product_name_parts['substituent_1']}-{product_name_parts['substituent_2']}{product_name_parts['parent_ring']}-{product_name_parts['unsaturation']}-{product_name_parts['principal_group']}"
# A cleaner representation of the name.
clean_name = "(E)-2-methoxy-10-methylcyclodec-2-en-1-one"

# Step 3: Print the final IUPAC name.
print(f"The IUPAC name of the product is: {clean_name}")

# Step 4: As requested, output each number from the name.
print("\nThe numbers in the IUPAC name are:")
numbers = [2, 10, 2, 1]
for number in numbers:
    print(number)
