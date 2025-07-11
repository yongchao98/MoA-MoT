# Number of atoms in toluene (C7H8)
num_carbon = 7
num_hydrogen = 8

# Number of contracted Gaussian functions for each atom in the 6-31G basis set
# For Carbon (C): 1 core (1s) + 4 inner valence (2s, 2p) + 4 outer valence (2s', 2p') = 9
functions_per_carbon = 9
# For Hydrogen (H): 1 inner valence (1s) + 1 outer valence (1s') = 2
functions_per_hydrogen = 2

# Calculate the total number of functions from each element
total_carbon_functions = num_carbon * functions_per_carbon
total_hydrogen_functions = num_hydrogen * functions_per_hydrogen

# Calculate the total number of contracted Gaussian functions for the molecule
total_functions = total_carbon_functions + total_hydrogen_functions

# Print the breakdown of the calculation and the final answer
print(f"Molecule: Toluene (C{num_carbon}H{num_hydrogen})")
print(f"Basis Set: 6-31G")
print("-" * 30)
print(f"Number of functions per Carbon atom: {functions_per_carbon}")
print(f"Number of functions per Hydrogen atom: {functions_per_hydrogen}")
print("-" * 30)
print(f"Total functions = (Number of C atoms * Functions per C) + (Number of H atoms * Functions per H)")
print(f"Total functions = ({num_carbon} * {functions_per_carbon}) + ({num_hydrogen} * {functions_per_hydrogen})")
print(f"Total functions = {total_carbon_functions} + {total_hydrogen_functions}")
print(f"Total functions = {total_functions}")
