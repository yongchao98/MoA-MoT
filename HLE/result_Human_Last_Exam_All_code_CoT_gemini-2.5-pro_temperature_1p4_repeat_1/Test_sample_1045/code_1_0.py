# Number of atoms in toluene (C7H8)
num_carbon = 7
num_hydrogen = 8

# Number of contracted Gaussian functions per atom for the 6-31G basis set
# For Carbon (C): 1 (1s) + 2 (2s) + 3*2 (2p) = 9
functions_per_carbon = 9
# For Hydrogen (H): 2 (1s)
functions_per_hydrogen = 2

# Calculate the total number of functions for each element
total_carbon_functions = num_carbon * functions_per_carbon
total_hydrogen_functions = num_hydrogen * functions_per_hydrogen

# Calculate the grand total
total_functions = total_carbon_functions + total_hydrogen_functions

# Print the breakdown of the calculation
print(f"For toluene (C{num_carbon}H{num_hydrogen}) with a 6-31G basis set:")
print(f"Number of functions for Carbon = {num_carbon} atoms * {functions_per_carbon} functions/atom = {total_carbon_functions}")
print(f"Number of functions for Hydrogen = {num_hydrogen} atoms * {functions_per_hydrogen} functions/atom = {total_hydrogen_functions}")
print(f"Total number of contracted Gaussian functions = {total_carbon_functions} + {total_hydrogen_functions} = {total_functions}")
