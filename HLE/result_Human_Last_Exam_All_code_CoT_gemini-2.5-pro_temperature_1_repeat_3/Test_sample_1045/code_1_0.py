# Number of atoms in toluene (C7H8)
num_carbon_atoms = 7
num_hydrogen_atoms = 8

# Number of contracted Gaussian functions per atom for the 6-31G basis set
# For Carbon (C): 1s (core) + 2s,2p (inner valence) + 2s,2p (outer valence) = 1 + 4 + 4 = 9
functions_per_carbon = 9
# For Hydrogen (H): 1s (inner valence) + 1s (outer valence) = 1 + 1 = 2
functions_per_hydrogen = 2

# Calculate the total number of functions
total_functions_carbon = num_carbon_atoms * functions_per_carbon
total_functions_hydrogen = num_hydrogen_atoms * functions_per_hydrogen
total_functions = total_functions_carbon + total_functions_hydrogen

# Print the final equation and the result
print(f"For toluene (C7H8) with a 6-31G basis set:")
print(f"Total contracted functions = (Number of C atoms * Functions per C) + (Number of H atoms * Functions per H)")
print(f"({num_carbon_atoms} * {functions_per_carbon}) + ({num_hydrogen_atoms} * {functions_per_hydrogen}) = {total_functions}")