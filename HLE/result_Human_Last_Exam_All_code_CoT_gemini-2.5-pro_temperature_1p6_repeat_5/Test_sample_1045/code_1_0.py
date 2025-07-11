# Number of atoms in toluene (C7H8)
num_carbon = 7
num_hydrogen = 8

# For a 6-31G basis set:
# Number of contracted functions for a Hydrogen atom (1s orbital split into 2)
functions_per_H = 2
# Number of contracted functions for a Carbon atom (1s core is 1 function, 2s and 2p valence are split into 2 functions each)
# 1 (core) + 4 (valence orbitals) * 2 = 9
functions_per_C = 9

# Calculate the total number of basis functions
total_functions = (num_carbon * functions_per_C) + (num_hydrogen * functions_per_H)

# Print the final equation and the result
print(f"For toluene (C7H8) with a 6-31G basis set, the calculation is:")
print(f"({num_carbon} C * {functions_per_C} functions/C) + ({num_hydrogen} H * {functions_per_H} functions/H) = {total_functions}")
print(f"Total number of contracted Gaussian functions = {total_functions}")