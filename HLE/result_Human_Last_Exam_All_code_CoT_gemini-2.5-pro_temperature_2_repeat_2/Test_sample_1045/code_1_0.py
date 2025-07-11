# Define the number of atoms in toluene (C7H8)
num_carbon_atoms = 7
num_hydrogen_atoms = 8

# Define the number of contracted Gaussian functions per atom for the 6-31G basis set
# For Hydrogen (1s orbital): split into 2 basis functions (3G, 1G)
functions_per_H = 2
# For Carbon (1s core, 2s, 2p valence): 1 core function (6G) + 2 sets of valence functions (4*(3G) + 4*(1G)) = 1 + 4 + 4 = 9
functions_per_C = 9

# Calculate the contribution from carbon and hydrogen atoms
carbon_contribution = num_carbon_atoms * functions_per_C
hydrogen_contribution = num_hydrogen_atoms * functions_per_H

# Calculate the total number of contracted Gaussian functions
total_functions = carbon_contribution + hydrogen_contribution

# Print the breakdown of the calculation and the final result
print(f"For a 6-31G basis set calculation of toluene (C{num_carbon_atoms}H{num_hydrogen_atoms}):")
print(f"Number of contracted functions per Carbon (C) atom: {functions_per_C}")
print(f"Number of contracted functions per Hydrogen (H) atom: {functions_per_H}")
print("\nFinal Calculation:")
print(f"Total contracted functions = ({num_carbon_atoms} * {functions_per_C}) + ({num_hydrogen_atoms} * {functions_per_H}) = {carbon_contribution} + {hydrogen_contribution} = {total_functions}")

<<<79>>>