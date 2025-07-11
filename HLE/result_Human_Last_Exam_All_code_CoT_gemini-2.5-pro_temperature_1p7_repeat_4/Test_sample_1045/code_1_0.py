# Number of atoms in toluene (C7H8)
num_C = 7
num_H = 8

# Number of contracted Gaussian functions per atom for the 6-31G basis set
# For Carbon (C): 1 (core 1s) + 2 (valence 2s) + 2*3 (valence 2p) = 9
funcs_per_C = 9
# For Hydrogen (H): 2 (split-valence 1s)
funcs_per_H = 2

# Calculate the total number of functions from each element
total_C = num_C * funcs_per_C
total_H = num_H * funcs_per_H

# Calculate the total for the entire molecule
total_funcs = total_C + total_H

# Print the breakdown of the calculation and the final answer
print(f"Molecule: Toluene (C{num_C}H{num_H})")
print(f"Basis Set: 6-31G")
print("-" * 30)
print(f"Functions per Carbon (C) atom: {funcs_per_C}")
print(f"Functions per Hydrogen (H) atom: {funcs_per_H}")
print("-" * 30)
print("The total number of contracted Gaussian functions is calculated as follows:")
print(f"Total = (Atoms of C * Functions per C) + (Atoms of H * Functions per H)")
# The final equation with each number and the final result
print(f"({num_C} * {funcs_per_C}) + ({num_H} * {funcs_per_H}) = {total_C} + {total_H} = {total_funcs}")
