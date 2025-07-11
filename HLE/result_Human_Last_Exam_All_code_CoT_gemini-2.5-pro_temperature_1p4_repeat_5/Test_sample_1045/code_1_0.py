# Toluene molecule: C7H8
num_carbon_atoms = 7
num_hydrogen_atoms = 8

# Number of contracted Gaussian functions for each atom with the 6-31G basis set
# For Carbon (C): 1s(core) + 2s(valence_inner) + 2p(valence_inner) + 2s(valence_outer) + 2p(valence_outer)
# This simplifies to 1 (for 1s) + 4 (for inner 2sp) + 4 (for outer 2sp) = 9 functions.
cgf_per_carbon = 9

# For Hydrogen (H): 1s(valence_inner) + 1s(valence_outer) = 2 functions.
cgf_per_hydrogen = 2

# Calculate the total contribution from Carbon atoms
total_cgf_carbon = num_carbon_atoms * cgf_per_carbon

# Calculate the total contribution from Hydrogen atoms
total_cgf_hydrogen = num_hydrogen_atoms * cgf_per_hydrogen

# Calculate the total number of functions for the molecule
total_cgf = total_cgf_carbon + total_cgf_hydrogen

# Print the breakdown of the calculation
print("Calculation for Toluene (C7H8) with 6-31G basis set:\n")
print(f"Number of Carbon atoms: {num_carbon_atoms}")
print(f"Number of Hydrogen atoms: {num_hydrogen_atoms}\n")
print(f"Functions per Carbon atom: {cgf_per_carbon}")
print(f"Functions per Hydrogen atom: {cgf_per_hydrogen}\n")

# Print the final equation with all numbers
print("Final Equation:")
print(f"Total functions = ({num_carbon_atoms} C atoms * {cgf_per_carbon} funcs/C) + ({num_hydrogen_atoms} H atoms * {cgf_per_hydrogen} funcs/H)")
print(f"                = {total_cgf_carbon} + {total_cgf_hydrogen}")
print(f"                = {total_cgf}")
