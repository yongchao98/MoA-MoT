import re

# This script provides the IUPAC name for the major product of the specified reaction
# and extracts the numerical subscripts from the balanced chemical equation.

# The overall balanced chemical equation for the two-step, one-pot reaction is:
# C13H18O2S -> C7H12O + C6H6SO
equation = "C13H18O2S -> C7H12O + C6H6SO"

# The major organic product (C7H12O) is a gamma,delta-unsaturated aldehyde.
# Its IUPAC name is derived from its structure: O=CH-CH2-CH2-CH=C(CH3)2.
iupac_name = "5,5-dimethylpent-4-enal"

print(f"The IUPAC name of the major product is: {iupac_name}")

# Using a regular expression to find all integer numbers in the equation string.
# These numbers correspond to the subscripts in the molecular formulas.
numbers_found = re.findall(r'\d+', equation)

# As requested, outputting each number from the final equation.
print(f"The numbers in the final balanced chemical equation ({equation}) are: {', '.join(numbers_found)}")