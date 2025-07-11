# Number of atoms in toluene (C7H8)
num_C = 7
num_H = 8

# Number of contracted Gaussian functions per atom for the 6-31G basis set
# For Carbon (C): 1 core (1s) + 4 valence-inner (2s, 2p) + 4 valence-outer (2s, 2p) = 9
functions_per_C = 9
# For Hydrogen (H): 1 valence-inner (1s) + 1 valence-outer (1s) = 2
functions_per_H = 2

# Calculate the total number of functions from Carbon atoms
total_functions_C = num_C * functions_per_C

# Calculate the total number of functions from Hydrogen atoms
total_functions_H = num_H * functions_per_H

# Calculate the total for the molecule
total_functions = total_functions_C + total_functions_H

# Print the breakdown and the final result
print("Calculation for Toluene (C7H8) with 6-31G basis set:")
print(f"Number of Carbon atoms: {num_C}")
print(f"Functions per Carbon atom: {functions_per_C}")
print(f"Number of Hydrogen atoms: {num_H}")
print(f"Functions per Hydrogen atom: {functions_per_H}")
print("\nTotal contracted Gaussian functions are calculated as:")
print(f"({num_C} C atoms * {functions_per_C} functions/C) + ({num_H} H atoms * {functions_per_H} functions/H)")
print(f"= {total_functions_C} + {total_functions_H}")
print(f"= {total_functions}")

print("\nFinal Equation:")
print(f"{num_C} * {functions_per_C} + {num_H} * {functions_per_H} = {total_functions}")