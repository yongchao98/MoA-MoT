# Number of atoms in toluene (C7H8)
num_carbon_atoms = 7
num_hydrogen_atoms = 8

# Number of contracted basis functions per atom for the 6-31G basis set
# For Carbon (C): 1 core (1s) + 2 valence (2s) + 3*2 valence (2p) = 9
functions_per_carbon = 9
# For Hydrogen (H): 2 valence (1s)
functions_per_hydrogen = 2

# Calculate contributions from each element
total_carbon_functions = num_carbon_atoms * functions_per_carbon
total_hydrogen_functions = num_hydrogen_atoms * functions_per_hydrogen

# Calculate the total number of basis functions
total_functions = total_carbon_functions + total_hydrogen_functions

# Print the detailed calculation and the final answer
print("Calculation for toluene (C7H8) with the 6-31G basis set:")
print(f"Number of Carbon atoms: {num_carbon_atoms}")
print(f"Number of Hydrogen atoms: {num_hydrogen_atoms}")
print(f"Basis functions per Carbon atom: {functions_per_carbon}")
print(f"Basis functions per Hydrogen atom: {functions_per_hydrogen}")
print("\nTotal contracted Gaussian functions = (C atoms * functions/C) + (H atoms * functions/H)")
print(f"Total = ({num_carbon_atoms} * {functions_per_carbon}) + ({num_hydrogen_atoms} * {functions_per_hydrogen})")
print(f"Total = {total_carbon_functions} + {total_hydrogen_functions}")
print(f"Total = {total_functions}")
