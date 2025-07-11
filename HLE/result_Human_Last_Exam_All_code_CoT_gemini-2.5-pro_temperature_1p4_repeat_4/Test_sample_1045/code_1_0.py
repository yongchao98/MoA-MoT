# Number of atoms in toluene (C7H8)
num_carbon = 7
num_hydrogen = 8

# Number of contracted Gaussian functions per atom for the 6-31G basis set
# For Carbon: 1 (core 1s) + 4 (inner valence 2s, 2p) + 4 (outer valence 2s, 2p) = 9
functions_per_carbon = 9
# For Hydrogen: 1 (inner valence 1s) + 1 (outer valence 1s) = 2
functions_per_hydrogen = 2

# Calculate the total number of functions from each element
total_carbon_functions = num_carbon * functions_per_carbon
total_hydrogen_functions = num_hydrogen * functions_per_hydrogen

# Calculate the grand total for the molecule
total_functions = total_carbon_functions + total_hydrogen_functions

# Print the breakdown and the final equation
print(f"Molecule: Toluene (C{num_carbon}H{num_hydrogen})")
print("Basis Set: 6-31G")
print("-" * 30)
print(f"Functions per Carbon atom: {functions_per_carbon}")
print(f"Functions per Hydrogen atom: {functions_per_hydrogen}")
print("-" * 30)
print("Total number of contracted Gaussian functions:")
# The final part of the code prints the equation as requested
print(f"({num_carbon} C atoms * {functions_per_carbon} functions/C) + ({num_hydrogen} H atoms * {functions_per_hydrogen} functions/H) = {total_carbon_functions} + {total_hydrogen_functions} = {total_functions}")
