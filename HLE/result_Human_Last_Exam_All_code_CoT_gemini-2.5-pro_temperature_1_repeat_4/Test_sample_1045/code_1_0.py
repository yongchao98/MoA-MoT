# Number of atoms in toluene (C7H8)
num_carbon_atoms = 7
num_hydrogen_atoms = 8

# Number of contracted Gaussian functions per atom for the 6-31G basis set
# For Carbon (C): 1 (1s) + 1 (2s_inner) + 3 (2p_inner) + 1 (2s_outer) + 3 (2p_outer) = 9
functions_per_carbon = 9
# For Hydrogen (H): 1 (1s_inner) + 1 (1s_outer) = 2
functions_per_hydrogen = 2

# Calculate the contribution from Carbon atoms
total_functions_carbon = num_carbon_atoms * functions_per_carbon

# Calculate the contribution from Hydrogen atoms
total_functions_hydrogen = num_hydrogen_atoms * functions_per_hydrogen

# Calculate the total number of contracted Gaussian functions
total_functions = total_functions_carbon + total_functions_hydrogen

# Print the final equation and the result
print(f"({num_carbon_atoms} C atoms * {functions_per_carbon} functions/C) + ({num_hydrogen_atoms} H atoms * {functions_per_hydrogen} functions/H) = {total_functions_carbon} + {total_functions_hydrogen} = {total_functions}")
