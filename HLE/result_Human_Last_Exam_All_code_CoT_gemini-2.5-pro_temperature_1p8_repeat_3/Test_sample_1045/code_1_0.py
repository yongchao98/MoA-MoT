# Toluene has the chemical formula C7H8.
num_carbon_atoms = 7
num_hydrogen_atoms = 8

# In the 6-31G basis set:
# Each Carbon (C) atom is described by 9 contracted Gaussian functions.
# (1 for the core 1s orbital, 4 for the inner valence shell, 4 for the outer valence shell)
functions_per_carbon = 9

# Each Hydrogen (H) atom is described by 2 contracted Gaussian functions.
# (1 for the inner valence shell, 1 for the outer valence shell)
functions_per_hydrogen = 2

# Calculate the contribution from Carbon and Hydrogen atoms
carbon_contribution = num_carbon_atoms * functions_per_carbon
hydrogen_contribution = num_hydrogen_atoms * functions_per_hydrogen

# Calculate the total number of functions
total_functions = carbon_contribution + hydrogen_contribution

# Print the explanation and the final calculation
print(f"The calculation for the total number of contracted Gaussian functions for toluene (C7H8) with a 6-31G basis set is as follows:")
print(f"({num_carbon_atoms} C atoms * {functions_per_carbon} functions/C) + ({num_hydrogen_atoms} H atoms * {functions_per_hydrogen} functions/H)")
print(f"= ({carbon_contribution}) + ({hydrogen_contribution})")
print(f"= {total_functions}")