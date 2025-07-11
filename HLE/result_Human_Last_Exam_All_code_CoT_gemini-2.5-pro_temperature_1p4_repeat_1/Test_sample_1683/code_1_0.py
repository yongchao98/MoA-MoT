# This script identifies the final product of the multi-step synthesis.
# The analysis shows the transformation of (2-bromophenyl)methanol through
# several intermediates to the final compound.

# The final step is the Jones oxidation of 1,2-bis(hydroxymethyl)benzene (Compound 3).
# Jones reagent oxidizes primary alcohols to carboxylic acids.
# C6H4(CH2OH)2 -> C6H4(COOH)2

# Define the properties of the final product, Compound 4.
compound_name = "Phthalic acid"
common_name_alt = "Benzene-1,2-dicarboxylic acid"

# Molecular formula of Phthalic Acid is C8H6O4.
molecular_formula = "C8H6O4"

# Define the number of atoms for each element in the final molecule's "equation".
carbon_atoms = 8
hydrogen_atoms = 6
oxygen_atoms = 4

# Print the final answer and the atomic composition.
print(f"The final product, Compound 4, is {compound_name}.")
print(f"The molecular formula is {molecular_formula}.")
print("\nThe composition of the final molecule is:")
print(f"Number of Carbon atoms: {carbon_atoms}")
print(f"Number of Hydrogen atoms: {hydrogen_atoms}")
print(f"Number of Oxygen atoms: {oxygen_atoms}")
