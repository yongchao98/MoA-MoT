# The molecule is toluene (C7H8)
num_carbon_atoms = 7
num_hydrogen_atoms = 8

# For the 6-31G basis set:
# Each Carbon atom contributes 9 contracted Gaussian functions.
# (1 for core 1s, and 2*4=8 for the split valence 2s and 2p orbitals)
funcs_per_carbon = 9

# Each Hydrogen atom contributes 2 contracted Gaussian functions.
# (The 1s orbital is split into two functions)
funcs_per_hydrogen = 2

# Calculate the total number of functions from carbon atoms
total_carbon_funcs = num_carbon_atoms * funcs_per_carbon

# Calculate the total number of functions from hydrogen atoms
total_hydrogen_funcs = num_hydrogen_atoms * funcs_per_hydrogen

# Calculate the grand total for the molecule
total_funcs = total_carbon_funcs + total_hydrogen_funcs

# Print the breakdown of the calculation and the final answer
print(f"To find the total number of contracted Gaussian functions for toluene (C7H8) with a 6-31G basis set, we calculate:")
print(f"({num_carbon_atoms} * {funcs_per_carbon}) + ({num_hydrogen_atoms} * {funcs_per_hydrogen}) = {total_funcs}")