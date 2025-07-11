# Number of atoms in toluene (C7H8)
num_C = 7
num_H = 8

# Number of contracted Gaussian functions for each atom with the 6-31G basis set
# For Carbon (C):
# 1 core function (1s)
# 4 inner valence functions (2s, 2px, 2py, 2pz)
# 4 outer valence functions (2s, 2px, 2py, 2pz)
# Total = 1 + 4 + 4 = 9
functions_per_C = 9

# For Hydrogen (H):
# No core electrons
# 1 inner valence function (1s)
# 1 outer valence function (1s)
# Total = 1 + 1 = 2
functions_per_H = 2

# Calculate the total number of functions from Carbon atoms
total_C_functions = num_C * functions_per_C

# Calculate the total number of functions from Hydrogen atoms
total_H_functions = num_H * functions_per_H

# Calculate the total for the molecule
total_functions = total_C_functions + total_H_functions

# Print the breakdown of the calculation and the final result
print(f"Molecule: Toluene (C{num_C}H{num_H})")
print(f"Basis Set: 6-31G")
print("-" * 30)
print(f"Functions per Carbon (C) atom: {functions_per_C}")
print(f"Functions per Hydrogen (H) atom: {functions_per_H}")
print("-" * 30)
print("Total number of contracted Gaussian functions:")
print(f"({num_C} * {functions_per_C}) + ({num_H} * {functions_per_H}) = {total_functions}")
