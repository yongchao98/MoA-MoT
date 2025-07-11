# Toluene chemical formula: C7H8
num_carbon_atoms = 7
num_hydrogen_atoms = 8

# Number of contracted Gaussian functions for each atom in the 6-31G basis set
# For Carbon (heavy atom): 1 core (1s) + 4 inner valence (2s, 2p) + 4 outer valence (2s, 2p) = 9
functions_per_carbon = 9
# For Hydrogen: 1 inner valence (1s) + 1 outer valence (1s) = 2
functions_per_hydrogen = 2

# Calculate the total number of functions for each element
total_functions_carbon = num_carbon_atoms * functions_per_carbon
total_functions_hydrogen = num_hydrogen_atoms * functions_per_hydrogen

# Calculate the total number of contracted Gaussian functions for the molecule
total_functions = total_functions_carbon + total_functions_hydrogen

# Print the breakdown and the final result
print(f"For a 6-31G basis set calculation of toluene (C{num_carbon_atoms}H{num_hydrogen_atoms}):")
print(f"Each of the {num_carbon_atoms} Carbon atoms uses {functions_per_carbon} contracted functions.")
print(f"Each of the {num_hydrogen_atoms} Hydrogen atoms uses {functions_per_hydrogen} contracted functions.")
print("\nThe total number of contracted Gaussian functions is calculated as:")
print(f"Total = ({num_carbon_atoms} Carbon atoms * {functions_per_carbon} functions/atom) + ({num_hydrogen_atoms} Hydrogen atoms * {functions_per_hydrogen} functions/atom)")
print(f"Total = ({total_functions_carbon}) + ({total_functions_hydrogen}) = {total_functions}")
