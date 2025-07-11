# Number of atoms in toluene (C7H8)
num_carbon_atoms = 7
num_hydrogen_atoms = 8

# Number of contracted Gaussian functions per atom for the 6-31G basis set
# For Carbon: 1 (core 1s) + 4 (inner valence 2s, 2p) + 4 (outer valence 2s, 2p) = 9
functions_per_carbon = 9
# For Hydrogen: 1 (inner valence 1s) + 1 (outer valence 1s) = 2
functions_per_hydrogen = 2

# Calculate the total number of contracted functions
total_functions = (num_carbon_atoms * functions_per_carbon) + (num_hydrogen_atoms * functions_per_hydrogen)

# Print the breakdown of the calculation and the final result
print(f"Molecule: Toluene (C{num_carbon_atoms}H{num_hydrogen_atoms})")
print("Basis Set: 6-31G")
print(f"Functions per Carbon atom: {functions_per_carbon}")
print(f"Functions per Hydrogen atom: {functions_per_hydrogen}")
print("\nTotal number of contracted functions calculation:")
print(f"({num_carbon_atoms} C atoms * {functions_per_carbon} functions/C) + ({num_hydrogen_atoms} H atoms * {functions_per_hydrogen} functions/H)")
print(f"= {num_carbon_atoms * functions_per_carbon} + {num_hydrogen_atoms * functions_per_hydrogen}")
print(f"= {total_functions}")
