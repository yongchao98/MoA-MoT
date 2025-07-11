# Number of atoms in toluene (C7H8)
num_carbon = 7
num_hydrogen = 8

# Number of contracted Gaussian functions per atom for the 6-31G basis set
# For Carbon: 1 (1s core) + 2 (2s valence) + 6 (2p valence) = 9
functions_per_carbon = 9
# For Hydrogen: 2 (1s valence)
functions_per_hydrogen = 2

# Calculate the total number of functions
total_functions_carbon = num_carbon * functions_per_carbon
total_functions_hydrogen = num_hydrogen * functions_per_hydrogen
total_functions = total_functions_carbon + total_functions_hydrogen

# Print the explanation and the final calculation
print("Calculation for the number of contracted Gaussian functions for toluene (C7H8) with the 6-31G basis set.")
print("-" * 80)
print(f"Number of Carbon atoms: {num_carbon}")
print(f"Number of Hydrogen atoms: {num_hydrogen}")
print("\nAccording to the 6-31G basis set:")
print(f"  - Each Carbon atom contributes {functions_per_carbon} functions.")
print(f"  - Each Hydrogen atom contributes {functions_per_hydrogen} functions.")
print("\nThe total number of functions is calculated as:")
print(f"Total = (Number of C atoms * Functions per C) + (Number of H atoms * Functions per H)")
print(f"Total = ({num_carbon} * {functions_per_carbon}) + ({num_hydrogen} * {functions_per_hydrogen})")
print(f"Total = {total_functions_carbon} + {total_functions_hydrogen}")
print(f"Total = {total_functions}")

print(f"\nFinal Answer: {total_functions} contracted Gaussian functions are used.")
<<<79>>>