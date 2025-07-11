# Number of atoms in toluene (C7H8)
num_carbon = 7
num_hydrogen = 8

# Number of contracted Gaussian functions per atom for the 6-31G basis set
# For Carbon (C): 1 (core 1s) + 4 (inner valence 2s, 2p) + 4 (outer valence 2s, 2p) = 9
funcs_per_carbon = 9
# For Hydrogen (H): The 1s orbital is split into two functions.
funcs_per_hydrogen = 2

# Calculate the total number of functions
total_funcs_carbon = num_carbon * funcs_per_carbon
total_funcs_hydrogen = num_hydrogen * funcs_per_hydrogen
total_functions = total_funcs_carbon + total_funcs_hydrogen

# Print the explanation and the final equation
print("For a 6-31G basis set calculation of toluene (C7H8):")
print(f"Each of the {num_carbon} Carbon atoms uses {funcs_per_carbon} contracted Gaussian functions.")
print(f"Each of the {num_hydrogen} Hydrogen atoms uses {funcs_per_hydrogen} contracted Gaussian functions.")
print("\nThe total number of contracted functions is calculated as:")
print(f"({num_carbon} C atoms * {funcs_per_carbon} functions/C) + ({num_hydrogen} H atoms * {funcs_per_hydrogen} functions/H)")
print(f"= {num_carbon} * {funcs_per_carbon} + {num_hydrogen} * {funcs_per_hydrogen} = {total_functions}")
