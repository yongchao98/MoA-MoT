# Number of atoms in toluene (C7H8)
num_carbon = 7
num_hydrogen = 8

# Number of contracted Gaussian functions per atom for the 6-31G basis set
# For Carbon: 1 (core) + 4 (inner valence s,p) + 4 (outer valence s,p) = 9
functions_per_carbon = 9
# For Hydrogen: 1 (inner valence s) + 1 (outer valence s) = 2
functions_per_hydrogen = 2

# Calculate the total number of functions from each element
total_carbon_functions = num_carbon * functions_per_carbon
total_hydrogen_functions = num_hydrogen * functions_per_hydrogen

# Calculate the grand total
total_functions = total_carbon_functions + total_hydrogen_functions

# Print the final calculation and result
print(f"For toluene (C7H8) with a 6-31G basis set:")
print(f"The total number of contracted basis functions is calculated as:")
print(f"({num_carbon} C atoms * {functions_per_carbon} functions/C) + ({num_hydrogen} H atoms * {functions_per_hydrogen} functions/H)")
print(f"= ({total_carbon_functions}) + ({total_hydrogen_functions})")
print(f"= {total_functions}")