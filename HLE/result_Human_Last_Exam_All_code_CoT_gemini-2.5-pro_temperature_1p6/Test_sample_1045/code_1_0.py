# This script calculates the total number of contracted Gaussian functions
# for a toluene molecule (C7H8) using the 6-31G basis set.

# Step 1: Define the number of atoms in toluene.
num_carbon_atoms = 7
num_hydrogen_atoms = 8

# Step 2: Define the number of contracted functions per atom for the 6-31G basis set.
# For Carbon (C): 1 (core) + 4 (inner valence) + 4 (outer valence) = 9 functions.
functions_per_carbon = 9
# For Hydrogen (H): 1 (inner valence) + 1 (outer valence) = 2 functions.
functions_per_hydrogen = 2

# Step 3: Calculate the total number of functions from each element.
total_carbon_functions = num_carbon_atoms * functions_per_carbon
total_hydrogen_functions = num_hydrogen_atoms * functions_per_hydrogen

# Step 4: Sum the functions to get the total for the molecule.
total_functions = total_carbon_functions + total_hydrogen_functions

# Step 5: Print the breakdown of the calculation and the final result.
print(f"To find the total number of contracted Gaussian functions for toluene (C{num_carbon_atoms}H{num_hydrogen_atoms}) with a 6-31G basis set:")
print(f"First, calculate the functions for the {num_carbon_atoms} Carbon atoms: {num_carbon_atoms} * {functions_per_carbon} = {total_carbon_functions}")
print(f"Next, calculate the functions for the {num_hydrogen_atoms} Hydrogen atoms: {num_hydrogen_atoms} * {functions_per_hydrogen} = {total_hydrogen_functions}")
print(f"Finally, sum them together for the total number of functions.")
print(f"Final Equation: {total_carbon_functions} + {total_hydrogen_functions} = {total_functions}")
