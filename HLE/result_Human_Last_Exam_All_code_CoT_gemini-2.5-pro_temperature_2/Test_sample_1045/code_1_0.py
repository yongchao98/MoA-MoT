# Number of atoms in Toluene (C7H8)
num_carbon_atoms = 7
num_hydrogen_atoms = 8

# Number of contracted Gaussian functions per atom for the 6-31G basis set
# For Carbon (C): 1 core (1s) + 4 inner valence (2s, 2p) + 4 outer valence (2s, 2p) = 9
functions_per_carbon = 9
# For Hydrogen (H): 1 inner valence (1s) + 1 outer valence (1s) = 2
functions_per_hydrogen = 2

# Calculate the total number of functions contributed by each element
total_from_carbon = num_carbon_atoms * functions_per_carbon
total_from_hydrogen = num_hydrogen_atoms * functions_per_hydrogen

# Calculate the grand total
total_functions = total_from_carbon + total_from_hydrogen

# Print the step-by-step calculation and the final result
print(f"For toluene (C{num_carbon_atoms}H{num_hydrogen_atoms}) with the 6-31G basis set:")
print(f"Each carbon atom contributes {functions_per_carbon} functions.")
print(f"Each hydrogen atom contributes {functions_per_hydrogen} functions.")
print("\nThe total number of contracted Gaussian functions is calculated as follows:")
print(f"({num_carbon_atoms} C atoms * {functions_per_carbon} functions/C) + ({num_hydrogen_atoms} H atoms * {functions_per_hydrogen} functions/H) = {total_functions}")
print(f"{total_from_carbon} + {total_from_hydrogen} = {total_functions}")
