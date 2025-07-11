# Number of atoms in toluene (C7H8)
num_carbon = 7
num_hydrogen = 8

# Number of contracted Gaussian functions for each atom in the 6-31G basis set
functions_per_carbon = 9  # 1 (1s) + 4 (2s, 2p) + 4 (2s', 2p')
functions_per_hydrogen = 2 # 1 (1s) + 1 (1s')

# Calculate the contribution from each type of atom
total_for_carbon = num_carbon * functions_per_carbon
total_for_hydrogen = num_hydrogen * functions_per_hydrogen

# Calculate the total number of basis functions
total_functions = total_for_carbon + total_for_hydrogen

# Print the breakdown of the calculation and the final result
print(f"Molecule: Toluene (C{num_carbon}H{num_hydrogen})")
print(f"Basis Set: 6-31G\n")
print(f"Functions per Carbon atom = {functions_per_carbon}")
print(f"Functions per Hydrogen atom = {functions_per_hydrogen}\n")
print(f"Contribution from {num_carbon} Carbon atoms: {num_carbon} * {functions_per_carbon} = {total_for_carbon}")
print(f"Contribution from {num_hydrogen} Hydrogen atoms: {num_hydrogen} * {functions_per_hydrogen} = {total_for_hydrogen}\n")
print(f"Total number of contracted Gaussian functions = ({num_carbon} * {functions_per_carbon}) + ({num_hydrogen} * {functions_per_hydrogen}) = {total_functions}")