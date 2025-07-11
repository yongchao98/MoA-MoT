# Number of atoms in toluene (C7H8)
num_carbon = 7
num_hydrogen = 8

# Number of contracted Gaussian functions per atom for the 6-31G basis set
# For Carbon (C): 1 core (1s) + 4 inner valence (2s, 2p) + 4 outer valence (2s, 2p) = 9
functions_per_carbon = 9
# For Hydrogen (H): 1 inner valence (1s) + 1 outer valence (1s) = 2
functions_per_hydrogen = 2

# Calculate the total number of functions from carbon atoms
total_from_carbon = num_carbon * functions_per_carbon

# Calculate the total number of functions from hydrogen atoms
total_from_hydrogen = num_hydrogen * functions_per_hydrogen

# Calculate the grand total
total_functions = total_from_carbon + total_from_hydrogen

# Print the breakdown of the calculation and the final result
print(f"For toluene (C7H8) with a 6-31G basis set:")
print(f"Number of Carbon atoms: {num_carbon}")
print(f"Number of Hydrogen atoms: {num_hydrogen}")
print(f"Functions per Carbon atom: {functions_per_carbon}")
print(f"Functions per Hydrogen atom: {functions_per_hydrogen}")
print("\nTotal contracted functions calculation:")
print(f"({num_carbon} * {functions_per_carbon}) + ({num_hydrogen} * {functions_per_hydrogen}) = {total_functions}")