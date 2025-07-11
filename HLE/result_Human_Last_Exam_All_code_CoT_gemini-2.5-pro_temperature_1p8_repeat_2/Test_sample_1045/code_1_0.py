# Number of atoms in toluene (C7H8)
num_carbon = 7
num_hydrogen = 8

# Number of contracted Gaussian functions per atom for the 6-31G basis set
functions_per_carbon = 9  # 1 (core) + 4 (inner valence) + 4 (outer valence)
functions_per_hydrogen = 2 # 1 (inner valence) + 1 (outer valence)

# Calculate the total number of contracted functions
total_functions = (num_carbon * functions_per_carbon) + (num_hydrogen * functions_per_hydrogen)

# Print the final equation and the result
print(f"The calculation for the total number of contracted Gaussian functions is:")
print(f"({num_carbon} Carbon atoms * {functions_per_carbon} functions/atom) + ({num_hydrogen} Hydrogen atoms * {functions_per_hydrogen} functions/atom)")
print(f"= {num_carbon * functions_per_carbon} + {num_hydrogen * functions_per_hydrogen}")
print(f"= {total_functions}")

print(f"\nFinal Equation:")
print(f"{num_carbon} * {functions_per_carbon} + {num_hydrogen} * {functions_per_hydrogen} = {total_functions}")