# Number of atoms in toluene (C7H8)
num_C = 7
num_H = 8

# Number of contracted Gaussian functions per atom for the 6-31G basis set
# For Carbon (C): 1 core (1s) + 4 inner valence (2s, 2p) + 4 outer valence (2s, 2p) = 9
functions_per_C = 9
# For Hydrogen (H): 1 inner valence (1s) + 1 outer valence (1s) = 2
functions_per_H = 2

# Calculate the total number of functions from Carbon atoms
total_C_functions = num_C * functions_per_C

# Calculate the total number of functions from Hydrogen atoms
total_H_functions = num_H * functions_per_H

# Calculate the grand total for the molecule
total_functions = total_C_functions + total_H_functions

# Print the breakdown of the calculation and the final result
print(f"Toluene (C7H8) with a 6-31G basis set:")
print(f"Number of Carbon atoms: {num_C}")
print(f"Number of Hydrogen atoms: {num_H}")
print("-" * 30)
print(f"Functions per Carbon atom: {functions_per_C}")
print(f"Functions per Hydrogen atom: {functions_per_H}")
print("-" * 30)
print(f"Total functions = (C atoms * functions/C) + (H atoms * functions/H)")
print(f"Total functions = ({num_C} * {functions_per_C}) + ({num_H} * {functions_per_H})")
print(f"Total functions = {total_C_functions} + {total_H_functions}")
print(f"Total number of contracted Gaussian functions = {total_functions}")
<<<79>>>