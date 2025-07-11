# Number of atoms in toluene (C7H8)
num_carbon_atoms = 7
num_hydrogen_atoms = 8

# Number of contracted Gaussian functions per atom for the 6-31G basis set
# For Carbon (C): 1 core (1s) + 4 inner valence (2s, 2p) + 4 outer valence (2s, 2p) = 9
functions_per_carbon = 9
# For Hydrogen (H): 1 inner valence (1s) + 1 outer valence (1s) = 2
functions_per_hydrogen = 2

# Calculate the total number of functions from carbon atoms
total_functions_carbon = num_carbon_atoms * functions_per_carbon

# Calculate the total number of functions from hydrogen atoms
total_functions_hydrogen = num_hydrogen_atoms * functions_per_hydrogen

# Calculate the grand total
total_functions = total_functions_carbon + total_functions_hydrogen

# Print the final calculation breakdown
print("Calculation for Toluene (C7H8) with 6-31G basis set:")
print(f"({num_carbon_atoms} C atoms * {functions_per_carbon} functions/C) + ({num_hydrogen_atoms} H atoms * {functions_per_hydrogen} functions/H) = {total_functions_carbon} + {total_functions_hydrogen} = {total_functions}")