# Number of atoms in toluene (C7H8)
num_C = 7
num_H = 8

# Number of contracted Gaussian functions per atom for the 6-31G basis set
funcs_per_C = 9  # 1 for core 1s, 4 for inner valence (2s,2p), 4 for outer valence (2s,2p)
funcs_per_H = 2  # 1 for inner valence 1s, 1 for outer valence 1s

# Calculate the total number of functions from each element
total_C_funcs = num_C * funcs_per_C
total_H_funcs = num_H * funcs_per_H

# Calculate the total for the molecule
total_funcs = total_C_funcs + total_H_funcs

print("To calculate the total number of contracted Gaussian functions for Toluene (C7H8) with the 6-31G basis set:")
print("We sum the functions for all atoms.")
print(f"Total = (Atoms of C * Functions per C) + (Atoms of H * Functions per H)")

# Print the final equation with the calculated numbers
print(f"Calculation: ({num_C} * {funcs_per_C}) + ({num_H} * {funcs_per_H}) = {total_funcs}")

print(f"\nThe total number of contracted Gaussian functions is: {total_funcs}")