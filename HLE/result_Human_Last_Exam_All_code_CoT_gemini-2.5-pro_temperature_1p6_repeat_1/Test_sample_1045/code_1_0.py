# Number of atoms in toluene (C7H8)
num_carbon_atoms = 7
num_hydrogen_atoms = 8

# Number of contracted Gaussian functions per atom for the 6-31G basis set
# For Carbon (C): 1 core (1s) + 4 inner valence (2s, 2p) + 4 outer valence (2s, 2p) = 9
funcs_per_carbon = 9
# For Hydrogen (H): 1 inner valence (1s) + 1 outer valence (1s) = 2
funcs_per_hydrogen = 2

# Calculate the contribution from each type of atom
carbon_contribution = num_carbon_atoms * funcs_per_carbon
hydrogen_contribution = num_hydrogen_atoms * funcs_per_hydrogen

# Calculate the total number of functions
total_functions = carbon_contribution + hydrogen_contribution

# Print the breakdown of the calculation
print(f"For a 6-31G basis set calculation of toluene (C{num_carbon_atoms}H{num_hydrogen_atoms}):")
print(f"Each of the {num_carbon_atoms} Carbon atoms uses {funcs_per_carbon} contracted functions.")
print(f"Each of the {num_hydrogen_atoms} Hydrogen atoms uses {funcs_per_hydrogen} contracted functions.")
print("\nCalculating the total:")
print(f"Total functions = ({num_carbon_atoms} C atoms * {funcs_per_carbon} funcs/C) + ({num_hydrogen_atoms} H atoms * {funcs_per_hydrogen} funcs/H)")
print(f"Total functions = {carbon_contribution} + {hydrogen_contribution}")
print(f"Total functions = {total_functions}")
