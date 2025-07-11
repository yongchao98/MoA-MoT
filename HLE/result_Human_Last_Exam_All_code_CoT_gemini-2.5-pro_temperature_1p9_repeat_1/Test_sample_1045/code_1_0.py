# Number of atoms in toluene (C7H8)
num_C_atoms = 7
num_H_atoms = 8

# Number of contracted Gaussian functions for each atom type in the 6-31G basis set
# For Carbon: 1 (core 1s) + 4 (inner valence 2s, 2p) + 4 (outer valence 2s, 2p) = 9
funcs_per_C = 9
# For Hydrogen: 1 (inner valence 1s) + 1 (outer valence 1s) = 2
funcs_per_H = 2

# Calculate the total number of functions contributed by each element
total_funcs_C = num_C_atoms * funcs_per_C
total_funcs_H = num_H_atoms * funcs_per_H

# Calculate the total number of functions for the molecule
total_funcs = total_funcs_C + total_funcs_H

# Print the final calculation and result
print(f"To calculate the total number of contracted Gaussian functions for toluene (C7H8) with a 6-31G basis set:")
print(f"Total = (Number of C atoms * Functions per C) + (Number of H atoms * Functions per H)")
print(f"Total = ({num_C_atoms} * {funcs_per_C}) + ({num_H_atoms} * {funcs_per_H})")
print(f"Total = {total_funcs_C} + {total_funcs_H}")
print(f"Total = {total_funcs}")
