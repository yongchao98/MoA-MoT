# Define the number of atoms for Toluene (C7H8)
num_C = 7
num_H = 8

# Define the number of contracted basis functions for each atom in the 6-31G basis set
# For Carbon (C): 1 core (1s) + 4 inner valence (2s, 2p) + 4 outer valence (2s, 2p) = 9
functions_per_C = 9

# For Hydrogen (H): 1 inner valence (1s) + 1 outer valence (1s) = 2
functions_per_H = 2

# Calculate the total number of functions from Carbon atoms
total_from_C = num_C * functions_per_C

# Calculate the total number of functions from Hydrogen atoms
total_from_H = num_H * functions_per_H

# Calculate the grand total for the molecule
total_functions = total_from_C + total_from_H

# Print the breakdown and the final equation
print(f"Toluene has {num_C} Carbon atoms and {num_H} Hydrogen atoms.")
print(f"In a 6-31G basis set:")
print(f"- Each Carbon atom contributes {functions_per_C} contracted Gaussian functions.")
print(f"- Each Hydrogen atom contributes {functions_per_H} contracted Gaussian functions.")
print("\nCalculating the total number of functions:")
print(f"({num_C} C atoms * {functions_per_C} functions/C) + ({num_H} H atoms * {functions_per_H} functions/H) = {total_functions}")
print(f"{total_from_C} + {total_from_H} = {total_functions}")
