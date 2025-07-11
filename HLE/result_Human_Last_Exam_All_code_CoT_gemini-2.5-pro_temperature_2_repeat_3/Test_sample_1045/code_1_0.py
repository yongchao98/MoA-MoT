# Toluene molecule is C7H8
# The basis set is 6-31G

# 1. Define the number of atoms
num_C_atoms = 7
num_H_atoms = 8

# 2. Define the number of contracted Gaussian functions per atom for the 6-31G basis set
# For a Carbon atom (C), it is 9 functions.
# For a Hydrogen atom (H), it is 2 functions.
funcs_per_C = 9
funcs_per_H = 2

# 3. Calculate the total number of functions
total_funcs_C = num_C_atoms * funcs_per_C
total_funcs_H = num_H_atoms * funcs_per_H
total_functions = total_funcs_C + total_funcs_H

# Print the final equation with all the numbers
print(f"The total number of contracted Gaussian functions for toluene (C7H8) with a 6-31G basis set is calculated as:")
print(f"({num_C_atoms} Carbon atoms * {funcs_per_C} functions) + ({num_H_atoms} Hydrogen atoms * {funcs_per_H} functions)")
print(f"Total = {num_C_atoms} * {funcs_per_C} + {num_H_atoms} * {funcs_per_H} = {total_functions}")
