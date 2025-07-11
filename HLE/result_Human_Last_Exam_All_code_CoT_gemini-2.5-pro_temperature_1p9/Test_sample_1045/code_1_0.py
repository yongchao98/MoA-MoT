# 1. Define the number of atoms in toluene (C7H8)
num_carbon_atoms = 7
num_hydrogen_atoms = 8

# 2. Define the number of contracted basis functions per atom for the 6-31G basis set
# For Carbon (C): 1 core (1s) + 4 inner valence (2s, 2p) + 4 outer valence (2s', 2p') = 9
functions_per_carbon = 9
# For Hydrogen (H): 1 inner valence (1s) + 1 outer valence (1s') = 2
functions_per_hydrogen = 2

# 3. Calculate the sub-totals for each element
total_carbon_functions = num_carbon_atoms * functions_per_carbon
total_hydrogen_functions = num_hydrogen_atoms * functions_per_hydrogen

# 4. Calculate the grand total
total_functions = total_carbon_functions + total_hydrogen_functions

# 5. Print the final calculation showing each step
print(f"To calculate the total contracted Gaussian functions for toluene (C7H8) with a 6-31G basis set:")
print(f"Total Functions = (Number of Carbon Atoms * Functions per C) + (Number of Hydrogen Atoms * Functions per H)")
print(f"Total Functions = ({num_carbon_atoms} * {functions_per_carbon}) + ({num_hydrogen_atoms} * {functions_per_hydrogen})")
print(f"Total Functions = {total_carbon_functions} + {total_hydrogen_functions}")
print(f"Total Functions = {total_functions}")