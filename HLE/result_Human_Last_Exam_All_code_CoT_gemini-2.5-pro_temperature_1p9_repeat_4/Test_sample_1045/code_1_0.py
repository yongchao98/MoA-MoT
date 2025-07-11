# Toluene molecule formula: C7H8
num_carbon_atoms = 7
num_hydrogen_atoms = 8

# Number of contracted Gaussian functions for each atom with the 6-31G basis set
# For Carbon (C): 1 core (1s) + 2 valence s (2s) + 2*3 valence p (2p) = 9
funcs_per_carbon = 9
# For Hydrogen (H): 2 valence s (1s) = 2
funcs_per_hydrogen = 2

# Calculate the total number of functions from carbon atoms
total_funcs_carbon = num_carbon_atoms * funcs_per_carbon

# Calculate the total number of functions from hydrogen atoms
total_funcs_hydrogen = num_hydrogen_atoms * funcs_per_hydrogen

# Calculate the total for the molecule
total_funcs_molecule = total_funcs_carbon + total_funcs_hydrogen

# Print the explanation and the final calculation
print("To calculate the total number of contracted Gaussian functions for Toluene (C7H8) with a 6-31G basis set:")
print(f"1. Count functions from the {num_carbon_atoms} Carbon atoms.")
print(f"   Each Carbon atom contributes {funcs_per_carbon} functions.")
print(f"   Subtotal for Carbon: {num_carbon_atoms} * {funcs_per_carbon} = {total_funcs_carbon}")
print("")
print(f"2. Count functions from the {num_hydrogen_atoms} Hydrogen atoms.")
print(f"   Each Hydrogen atom contributes {funcs_per_hydrogen} functions.")
print(f"   Subtotal for Hydrogen: {num_hydrogen_atoms} * {funcs_per_hydrogen} = {total_funcs_hydrogen}")
print("")
print("3. Sum the contributions from all atoms.")
print(f"   Total Functions = {total_funcs_carbon} (from C) + {total_funcs_hydrogen} (from H) = {total_funcs_molecule}")
print("\nFinal Answer:")
print(f"The total number of contracted Gaussian functions is {total_funcs_molecule}.")
