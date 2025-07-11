# Number of atoms in toluene (C7H8)
num_carbon_atoms = 7
num_hydrogen_atoms = 8

# Number of contracted Gaussian functions per atom for the 6-31G basis set
functions_per_carbon = 9  # 1s (core) + 2s (split) + 2p (split) -> 1 + 2 + (2*3) = 9
functions_per_hydrogen = 2 # 1s (split) -> 2

# Calculate the total number of functions for each element
total_functions_carbon = num_carbon_atoms * functions_per_carbon
total_functions_hydrogen = num_hydrogen_atoms * functions_per_hydrogen

# Calculate the grand total for the molecule
total_functions = total_functions_carbon + total_functions_hydrogen

# Print the breakdown of the calculation and the final result
print("To calculate the total number of contracted Gaussian functions for toluene (C7H8) with a 6-31G basis set:")
print(f"1. Carbon atoms: {num_carbon_atoms} atoms * {functions_per_carbon} functions/atom = {total_functions_carbon} functions")
print(f"2. Hydrogen atoms: {num_hydrogen_atoms} atoms * {functions_per_hydrogen} functions/atom = {total_functions_hydrogen} functions")
print("\nFinal Equation:")
print(f"({num_carbon_atoms} * {functions_per_carbon}) + ({num_hydrogen_atoms} * {functions_per_hydrogen}) = {total_functions_carbon} + {total_functions_hydrogen} = {total_functions}")
print(f"\nTotal contracted Gaussian functions = {total_functions}")