# Number of atoms in toluene (C7H8)
num_carbon_atoms = 7
num_hydrogen_atoms = 8

# Number of contracted Gaussian functions per atom for the 6-31G basis set
# For Carbon (a second-row atom): 1 for core (1s) + 4 for inner valence (2s, 2p) + 4 for outer valence (2s, 2p) = 9
funcs_per_carbon = 9
# For Hydrogen (a first-row atom): 1 for inner valence (1s) + 1 for outer valence (1s) = 2
funcs_per_hydrogen = 2

# Calculate the total number of functions from carbon and hydrogen atoms
total_funcs_carbon = num_carbon_atoms * funcs_per_carbon
total_funcs_hydrogen = num_hydrogen_atoms * funcs_per_hydrogen

# Calculate the grand total
total_contracted_functions = total_funcs_carbon + total_funcs_hydrogen

# Print the breakdown of the calculation and the final answer
print("Calculation for toluene (C7H8) with the 6-31G basis set:")
print(f"Total functions = ({num_carbon_atoms} C atoms * {funcs_per_carbon} functions/C) + ({num_hydrogen_atoms} H atoms * {funcs_per_hydrogen} functions/H)")
print(f"Total functions = {total_funcs_carbon} + {total_funcs_hydrogen}")
print(f"Total number of contracted Gaussian functions = {total_contracted_functions}")