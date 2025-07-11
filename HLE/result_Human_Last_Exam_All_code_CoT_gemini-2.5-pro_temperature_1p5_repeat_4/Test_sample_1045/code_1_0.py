# Define the number of atoms in toluene (C7H8)
num_carbon_atoms = 7
num_hydrogen_atoms = 8

# Define the number of contracted Gaussian functions per atom for the 6-31G basis set
# For Carbon: 1 (1s core) + 2 (2s valence) + 6 (2p valence) = 9
funcs_per_carbon = 9
# For Hydrogen: 2 (1s valence)
funcs_per_hydrogen = 2

# Calculate the total number of functions for each element
total_funcs_carbon = num_carbon_atoms * funcs_per_carbon
total_funcs_hydrogen = num_hydrogen_atoms * funcs_per_hydrogen

# Calculate the total for the molecule
total_contracted_functions = total_funcs_carbon + total_funcs_hydrogen

# Print the final result showing the equation
print(f"To calculate the total number of contracted Gaussian functions for toluene (C7H8) with the 6-31G basis set:")
print(f"First, calculate the contribution from the {num_carbon_atoms} Carbon atoms:")
print(f"    {num_carbon_atoms} C atoms * {funcs_per_carbon} functions/atom = {total_funcs_carbon} functions")
print(f"Next, calculate the contribution from the {num_hydrogen_atoms} Hydrogen atoms:")
print(f"    {num_hydrogen_atoms} H atoms * {funcs_per_hydrogen} functions/atom = {total_funcs_hydrogen} functions")
print("\nFinally, sum the contributions to get the total:")
print(f"Total Functions = ({num_carbon_atoms} * {funcs_per_carbon}) + ({num_hydrogen_atoms} * {funcs_per_hydrogen}) = {total_funcs_carbon} + {total_funcs_hydrogen} = {total_contracted_functions}")