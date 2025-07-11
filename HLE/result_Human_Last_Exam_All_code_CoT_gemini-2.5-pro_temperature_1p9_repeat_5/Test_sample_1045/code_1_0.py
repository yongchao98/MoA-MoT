# Define the number of atoms in toluene (C7H8)
num_carbon_atoms = 7
num_hydrogen_atoms = 8

# Define the number of contracted Gaussian functions per atom for the 6-31G basis set
# For Carbon: 1 core function + 4 inner valence functions + 4 outer valence functions = 9
funcs_per_carbon = 9

# For Hydrogen: 1 inner valence function + 1 outer valence function = 2
funcs_per_hydrogen = 2

# Calculate the total number of functions from Carbon atoms
total_funcs_from_carbon = num_carbon_atoms * funcs_per_carbon

# Calculate the total number of functions from Hydrogen atoms
total_funcs_from_hydrogen = num_hydrogen_atoms * funcs_per_hydrogen

# Calculate the total number of contracted Gaussian functions for the molecule
total_functions = total_funcs_from_carbon + total_funcs_from_hydrogen

# Print the final equation with each number explicitly shown
print("The total number of contracted Gaussian functions is calculated as:")
print(f"({num_carbon_atoms} Carbon atoms * {funcs_per_carbon} functions) + ({num_hydrogen_atoms} Hydrogen atoms * {funcs_per_hydrogen} functions)")
print(f"This gives the final equation:")
print(f"{num_carbon_atoms} * {funcs_per_carbon} + {num_hydrogen_atoms} * {funcs_per_hydrogen} = {total_functions}")
