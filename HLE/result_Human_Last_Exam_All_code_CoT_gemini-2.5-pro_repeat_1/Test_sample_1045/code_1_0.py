# Number of atoms in toluene (C7H8)
num_carbon_atoms = 7
num_hydrogen_atoms = 8

# Number of contracted basis functions per atom for the 6-31G basis set
functions_per_carbon = 9  # 1s (core) + 2s' (valence) + 2p' (valence) + 2s'' (valence) + 2p'' (valence) -> 1 + 1 + 3 + 1 + 3 = 9
functions_per_hydrogen = 2 # 1s' (valence) + 1s'' (valence) -> 1 + 1 = 2

# Calculate the total number of functions from each element
total_functions_carbon = num_carbon_atoms * functions_per_carbon
total_functions_hydrogen = num_hydrogen_atoms * functions_per_hydrogen

# Calculate the total for the molecule
total_functions = total_functions_carbon + total_functions_hydrogen

# Print the breakdown of the calculation and the final result
print(f"For a 6-31G basis set calculation of toluene (C{num_carbon_atoms}H{num_hydrogen_atoms}):")
print(f"Number of functions for {num_carbon_atoms} Carbon atoms = {num_carbon_atoms} * {functions_per_carbon} = {total_functions_carbon}")
print(f"Number of functions for {num_hydrogen_atoms} Hydrogen atoms = {num_hydrogen_atoms} * {functions_per_hydrogen} = {total_functions_hydrogen}")
print(f"Total number of contracted Gaussian functions = {total_functions_carbon} + {total_functions_hydrogen} = {total_functions}")
