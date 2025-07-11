# Number of atoms for toluene (C7H8)
num_carbon_atoms = 7
num_hydrogen_atoms = 8

# Number of contracted Gaussian functions per atom for the 6-31G basis set.
# For Carbon: 1 (core) + 4 (inner valence) + 4 (outer valence) = 9
funcs_per_carbon = 9
# For Hydrogen: 1 (inner valence) + 1 (outer valence) = 2
funcs_per_hydrogen = 2

# Calculate the total contributions from Carbon and Hydrogen atoms
total_funcs_from_carbon = num_carbon_atoms * funcs_per_carbon
total_funcs_from_hydrogen = num_hydrogen_atoms * funcs_per_hydrogen

# Calculate the total number of contracted Gaussian functions for the molecule
total_functions = total_funcs_from_carbon + total_funcs_from_hydrogen

# Print the final equation with all the numbers
print("The calculation for the total number of contracted Gaussian functions is:")
print(f"({num_carbon_atoms} Carbon atoms * {funcs_per_carbon} functions/atom) + ({num_hydrogen_atoms} Hydrogen atoms * {funcs_per_hydrogen} functions/atom)")
print(f"= {total_funcs_from_carbon} + {total_funcs_from_hydrogen}")
print(f"= {total_functions}")

print(f"\nFor a 6-31G basis set calculation of toluene (C7H8), a total of {total_functions} contracted Gaussian functions are used.")
