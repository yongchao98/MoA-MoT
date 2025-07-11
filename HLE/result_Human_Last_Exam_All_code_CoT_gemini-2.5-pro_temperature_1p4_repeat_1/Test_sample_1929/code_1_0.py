import re

# This script identifies the product of the Diels-Alder reaction between
# butadiene and 1,1-dichloro-2,2-difluoroethene and presents the reaction equation.

# Step 1: Define the IUPAC names of the reactants and the resulting product.
diene = "Butadiene"
dienophile = "1,1-dichloro-2,2-difluoroethene"
product = "3,3-dichloro-4,4-difluorocyclohex-1-ene"

# Step 2: Formulate the reaction equation as a string.
# Butadiene + 1,1-dichloro-2,2-difluoroethene -> 3,3-dichloro-4,4-difluorocyclohex-1-ene
# The stoichiometric coefficients are all 1.
final_equation = f"{diene} + {dienophile} -> {product}"

# Step 3: Print the full chemical equation.
print("The reaction between butadiene and 1,1-dichloro-2,2-difluoroethene is a Diels-Alder reaction.")
print("\nThe final equation describing this reaction is:")
print(final_equation)

# Step 4: As requested, extract and output each number present in the final equation.
# The numbers are part of the IUPAC names for the dienophile and the product.
numbers_in_equation = re.findall(r'\d', final_equation)

print("\nThe numbers from the chemical names in the final equation are:")
for number in numbers_in_equation:
    print(number)