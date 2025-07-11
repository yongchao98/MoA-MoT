# Toluene molecule information
num_carbon_atoms = 7
num_hydrogen_atoms = 8

# Number of contracted basis functions per atom for the 6-31G basis set
# For Hydrogen (1s valence orbital split in two): 2 functions
functions_per_hydrogen = 2

# For Carbon (1s core orbital + 2 sets of 2s, 2p valence orbitals):
# Core: 1 function (1s)
# Valence: 2 * (1 s-function + 3 p-functions) = 8 functions
# Total: 1 + 8 = 9 functions
functions_per_carbon = 9

# Calculate the total number of functions for each element
total_carbon_functions = num_carbon_atoms * functions_per_carbon
total_hydrogen_functions = num_hydrogen_atoms * functions_per_hydrogen

# Calculate the total for the molecule
total_functions = total_carbon_functions + total_hydrogen_functions

# Print the step-by-step explanation and the final result
print(f"To calculate the number of contracted Gaussian functions for Toluene (C{num_carbon_atoms}H{num_hydrogen_atoms}) with a 6-31G basis set:")
print("-" * 40)
print(f"1. Basis functions per Carbon (C) atom = {functions_per_carbon}")
print(f"2. Basis functions per Hydrogen (H) atom = {functions_per_hydrogen}")
print("-" * 40)
print(f"Total from {num_carbon_atoms} Carbon atoms = {num_carbon_atoms} * {functions_per_carbon} = {total_carbon_functions}")
print(f"Total from {num_hydrogen_atoms} Hydrogen atoms = {num_hydrogen_atoms} * {functions_per_hydrogen} = {total_hydrogen_functions}")
print("-" * 40)
print(f"Total contracted functions = {total_carbon_functions} + {total_hydrogen_functions} = {total_functions}")
<<<79>>>