# Define the number of atoms in toluene (C7H8)
num_carbon_atoms = 7
num_hydrogen_atoms = 8

# Define the number of contracted basis functions per atom for the 6-31G basis set
# For a second-row atom like Carbon: 1 core function + 2 sets of valence (s+p) functions = 1 + (1+3) + (1+3) = 9
funcs_per_carbon = 9
# For a Hydrogen atom: 2 sets of valence (s) functions = 1 + 1 = 2
funcs_per_hydrogen = 2

# Calculate the total number of functions from each element
total_funcs_carbon = num_carbon_atoms * funcs_per_carbon
total_funcs_hydrogen = num_hydrogen_atoms * funcs_per_hydrogen

# Calculate the total number of contracted Gaussian functions for the molecule
total_functions = total_funcs_carbon + total_funcs_hydrogen

# Print the final breakdown and the equation
print(f"Molecule: Toluene (C{num_carbon_atoms}H{num_hydrogen_atoms})")
print("Basis Set: 6-31G")
print("-" * 30)
print(f"Functions per Carbon atom: {funcs_per_carbon}")
print(f"Functions per Hydrogen atom: {funcs_per_hydrogen}")
print("-" * 30)
print("Total number of contracted Gaussian functions is calculated as:")
print(f"({num_carbon_atoms} C atoms * {funcs_per_carbon} funcs/atom) + ({num_hydrogen_atoms} H atoms * {funcs_per_hydrogen} funcs/atom) = {total_functions}")
print(f"So, the final equation is: {total_funcs_carbon} + {total_funcs_hydrogen} = {total_functions}")