# Number of atoms in toluene (C7H8)
num_carbon = 7
num_hydrogen = 8

# Number of contracted Gaussian functions for each atom in a 6-31G basis set
# For Carbon (C): 1s(core) + 2s(valence-inner) + 2s(valence-outer) + 2p(valence-inner)x3 + 2p(valence-outer)x3 = 1+1+1+3+3 = 9
# A simpler way: 1 core function + 4 valence functions (s, px, py, pz) * 2 (for the split) = 1 + 8 = 9
functions_per_carbon = 9

# For Hydrogen (H): 1s(valence-inner) + 1s(valence-outer) = 1 + 1 = 2
functions_per_hydrogen = 2

# Calculate the total number of basis functions
total_functions = (num_carbon * functions_per_carbon) + (num_hydrogen * functions_per_hydrogen)

# Print the breakdown of the calculation and the final result
print("Calculation for Toluene (C7H8) with a 6-31G basis set:")
print(f"Number of Carbon atoms: {num_carbon}")
print(f"Number of Hydrogen atoms: {num_hydrogen}")
print(f"Basis functions per Carbon atom: {functions_per_carbon}")
print(f"Basis functions per Hydrogen atom: {functions_per_hydrogen}")
print("\nTotal contracted Gaussian functions = (C atoms * functions/C) + (H atoms * functions/H)")
print(f"Total = ({num_carbon} * {functions_per_carbon}) + ({num_hydrogen} * {functions_per_hydrogen}) = {total_functions}")
