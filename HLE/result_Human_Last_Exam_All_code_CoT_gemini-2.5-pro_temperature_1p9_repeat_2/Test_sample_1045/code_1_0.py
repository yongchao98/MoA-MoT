# Toluene molecule composition
num_carbon = 7
num_hydrogen = 8

# Number of contracted Gaussian functions for each atom type with a 6-31G basis set
# For Carbon (C): 1s(core) + 2s(valence) + 2p(valence) -> 1 + 2 + (2*3) = 9 functions
# For Hydrogen (H): 1s(valence) -> 2 functions
funcs_per_carbon = 9
funcs_per_hydrogen = 2

# Calculate the total number of functions from Carbon atoms
total_carbon_funcs = num_carbon * funcs_per_carbon

# Calculate the total number of functions from Hydrogen atoms
total_hydrogen_funcs = num_hydrogen * funcs_per_hydrogen

# Calculate the total for the molecule
total_funcs = total_carbon_funcs + total_hydrogen_funcs

# Print the breakdown of the calculation
print(f"For toluene (C{num_carbon}H{num_hydrogen}) with a 6-31G basis set:")
print(f"Number of basis functions per Carbon atom = {funcs_per_carbon}")
print(f"Number of basis functions per Hydrogen atom = {funcs_per_hydrogen}")
print("\nCalculating the total number of contracted Gaussian functions:")
print(f"Total Functions = ({num_carbon} C atoms * {funcs_per_carbon} funcs/C) + ({num_hydrogen} H atoms * {funcs_per_hydrogen} funcs/H)")
print(f"Total Functions = ({num_carbon} * {funcs_per_carbon}) + ({num_hydrogen} * {funcs_per_hydrogen})")
print(f"Total Functions = {total_carbon_funcs} + {total_hydrogen_funcs}")
print(f"Total Functions = {total_funcs}")
