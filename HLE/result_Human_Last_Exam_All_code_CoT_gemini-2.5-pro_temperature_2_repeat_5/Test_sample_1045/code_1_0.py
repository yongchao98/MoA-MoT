# Number of atoms in toluene (C7H8)
num_C = 7
num_H = 8

# Number of contracted basis functions per atom for the 6-31G basis set
functions_per_C = 9  # 1 for core 1s, 2 for valence 2s, 6 for valence 2p (3*2)
functions_per_H = 2  # 2 for valence 1s

# Calculate the total number of functions from each element
total_C_functions = num_C * functions_per_C
total_H_functions = num_H * functions_per_H

# Calculate the total for the molecule
total_functions = total_C_functions + total_H_functions

# Print the explanation and the final equation
print("Calculation of contracted Gaussian functions for Toluene (C7H8) with a 6-31G basis set:")
print(f"Number of Carbon atoms: {num_C}")
print(f"Number of Hydrogen atoms: {num_H}")
print(f"Basis functions per Carbon atom: {functions_per_C}")
print(f"Basis functions per Hydrogen atom: {functions_per_H}")
print("\nThe total number of contracted Gaussian functions is calculated as:")
print(f"({num_C} C atoms * {functions_per_C} functions/C) + ({num_H} H atoms * {functions_per_H} functions/H) = {total_functions}")
print(f"This breaks down to: {total_C_functions} + {total_H_functions} = {total_functions}")
