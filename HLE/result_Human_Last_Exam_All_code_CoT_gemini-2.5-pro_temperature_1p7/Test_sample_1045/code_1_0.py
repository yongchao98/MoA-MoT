# Toluene molecule formula: C7H8
num_carbon_atoms = 7
num_hydrogen_atoms = 8

# Determine the number of contracted Gaussian functions for a 6-31G basis set.

# For Carbon (a second-row atom):
# Core (1s orbital): Represented by 1 contracted function (from the '6' in 6-31G).
# Valence (2s, 2p orbitals): Are split-valence (from the '31' in 6-31G).
# This means the 2s orbital gets 2 functions, and the three 2p orbitals also get 2 functions each.
# Total for one Carbon = 1 (for 1s) + 2 (for 2s) + 3*2 (for 2p_x, 2p_y, 2p_z) = 9
functions_per_carbon = 9

# For Hydrogen (a first-row atom):
# It only has a valence 1s orbital, which is split-valence.
# This means the 1s orbital gets 2 functions.
# Total for one Hydrogen = 2
functions_per_hydrogen = 2

# Calculate the total number of functions for the whole molecule
total_carbon_functions = num_carbon_atoms * functions_per_carbon
total_hydrogen_functions = num_hydrogen_atoms * functions_per_hydrogen
total_functions = total_carbon_functions + total_hydrogen_functions

# Print the breakdown and the final equation
print(f"Molecule: Toluene (C{num_carbon_atoms}H{num_hydrogen_atoms})")
print("Basis Set: 6-31G\n")
print(f"Number of functions for each Carbon atom: {functions_per_carbon}")
print(f"Number of functions for each Hydrogen atom: {functions_per_hydrogen}\n")
print("The final calculation is:")
print(f"Total functions = ({num_carbon_atoms} Carbon atoms * {functions_per_carbon} functions/atom) + ({num_hydrogen_atoms} Hydrogen atoms * {functions_per_hydrogen} functions/atom)")
print(f"Total functions = ({total_carbon_functions}) + ({total_hydrogen_functions})")
print(f"Total functions = {total_functions}")