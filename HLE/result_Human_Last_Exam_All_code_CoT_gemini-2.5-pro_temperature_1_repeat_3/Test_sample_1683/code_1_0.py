import re

# The final product of the reaction sequence is Phthalic acid.
compound_name = "Phthalic acid"
molecular_formula = "C8H6O4"

print(f"The final product, Compound 4, is {compound_name}.")
print(f"Its molecular formula is {molecular_formula}.")
print("\nThe numbers in the final 'equation' (molecular formula) are:")

# Use regex to find all sequences of digits in the formula string.
numbers = re.findall(r'\d+', molecular_formula)

# Print each number found.
for number in numbers:
    print(number)