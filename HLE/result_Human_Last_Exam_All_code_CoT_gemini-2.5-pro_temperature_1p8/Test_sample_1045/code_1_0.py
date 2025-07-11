# Number of atoms in toluene (C7H8)
num_C = 7
num_H = 8

# Number of contracted Gaussian functions per atom for the 6-31G basis set
# For Hydrogen (H): The 1s orbital is split into 2 functions (3G, 1G).
funcs_per_H = 2

# For Carbon (C):
# Core (1s): Represented by 1 function (6G).
# Valence (2s, 2p): Each of the 4 valence orbitals is split into 2 functions (3G, 1G).
# Total for C = 1 (core) + 4 * 2 (valence) = 9 functions.
funcs_per_C = 9

# Calculate the total number of functions from Carbon atoms
total_C = num_C * funcs_per_C

# Calculate the total number of functions from Hydrogen atoms
total_H = num_H * funcs_per_H

# Calculate the total for the entire molecule
total_funcs = total_C + total_H

# Print the explanation and the final equation
print("Calculation for toluene (C7H8) with 6-31G basis set:")
print(f"Number of Carbon atoms: {num_C}")
print(f"Number of Hydrogen atoms: {num_H}")
print(f"Contracted functions per Carbon atom: {funcs_per_C}")
print(f"Contracted functions per Hydrogen atom: {funcs_per_H}")
print("\nTotal number of contracted Gaussian functions:")
print(f"({num_C} * {funcs_per_C}) + ({num_H} * {funcs_per_H}) = {total_C} + {total_H} = {total_funcs}")
