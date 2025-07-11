# Toluene molecule: C7H8
# Basis set: 6-31G

# 1. Define the number of atoms for each element
num_carbon_atoms = 7
num_hydrogen_atoms = 8

# 2. Define the number of contracted basis functions per atom for the 6-31G basis set
# For Carbon (C): 1s (core) + 2s (valence) + 2p (valence)
# 1s is represented by 1 function.
# 2s is represented by 2 functions (split-valence).
# Each of the three 2p orbitals is represented by 2 functions (split-valence).
# Total for C = 1 (1s) + 2 (2s) + 3 * 2 (2p) = 9 functions.
functions_per_carbon = 9

# For Hydrogen (H): 1s (valence)
# 1s is represented by 2 functions (split-valence).
# Total for H = 2 functions.
functions_per_hydrogen = 2

# 3. Calculate the total number of contracted functions for the molecule
total_functions_carbon = num_carbon_atoms * functions_per_carbon
total_functions_hydrogen = num_hydrogen_atoms * functions_per_hydrogen
total_functions = total_functions_carbon + total_functions_hydrogen

# 4. Print the calculation steps and the final result
print("Calculation for Toluene (C7H8) with 6-31G basis set:")
print(f"Number of Carbon atoms: {num_carbon_atoms}")
print(f"Number of Hydrogen atoms: {num_hydrogen_atoms}")
print("-" * 20)
print(f"Basis functions per Carbon atom: {functions_per_carbon}")
print(f"Basis functions per Hydrogen atom: {functions_per_hydrogen}")
print("-" * 20)
print("Total functions calculation:")
print(f"({num_carbon_atoms} C atoms * {functions_per_carbon} functions/C) + ({num_hydrogen_atoms} H atoms * {functions_per_hydrogen} functions/H)")
print(f"= {total_functions_carbon} + {total_functions_hydrogen}")
print(f"= {total_functions}")
print("\nFinal Answer:")
print(f"Total contracted Gaussian functions = {total_functions}")