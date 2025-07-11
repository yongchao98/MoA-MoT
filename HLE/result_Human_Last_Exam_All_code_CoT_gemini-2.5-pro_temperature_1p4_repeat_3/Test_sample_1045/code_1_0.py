# Number of atoms in toluene (C7H8)
num_carbon = 7
num_hydrogen = 8

# For a 6-31G basis set:
# Number of contracted Gaussian functions for a Carbon atom is 9.
# (1 for the 1s core orbital, and 2 for each of the four valence orbitals: 2s, 2px, 2py, 2pz)
functions_per_carbon = 9

# Number of contracted Gaussian functions for a Hydrogen atom is 2.
# (The 1s valence orbital is split into 2 functions)
functions_per_hydrogen = 2

# Calculate the total number of functions from each element
total_carbon_functions = num_carbon * functions_per_carbon
total_hydrogen_functions = num_hydrogen * functions_per_hydrogen

# Calculate the final total
total_functions = total_carbon_functions + total_hydrogen_functions

# Print the breakdown of the calculation and the final result
print("The total number of contracted Gaussian functions is calculated as follows:")
print(f"({num_carbon} C atoms * {functions_per_carbon} functions/atom) + ({num_hydrogen} H atoms * {functions_per_hydrogen} functions/atom) = {total_carbon_functions} + {total_hydrogen_functions} = {total_functions}")