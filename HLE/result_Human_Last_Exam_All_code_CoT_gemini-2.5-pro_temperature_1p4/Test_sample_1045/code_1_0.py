# Toluene (C7H8) atom counts
num_carbon = 7
num_hydrogen = 8

# For the 6-31G basis set, the number of contracted Gaussian functions (basis functions) per atom is:
# - For Carbon (a second-row element): 9
# - For Hydrogen (a first-row element): 2
functions_per_carbon = 9
functions_per_hydrogen = 2

# Calculate the total number of functions
total_functions = (num_carbon * functions_per_carbon) + (num_hydrogen * functions_per_hydrogen)

# Print the breakdown and the final calculation
print("The total number of contracted Gaussian functions for Toluene (C7H8) with the 6-31G basis set is calculated by:")
print(f"Summing the functions from {num_carbon} Carbon atoms and {num_hydrogen} Hydrogen atoms.")
print(f"\n- Functions from Carbon: {num_carbon} atoms * {functions_per_carbon} functions/atom")
print(f"- Functions from Hydrogen: {num_hydrogen} atoms * {functions_per_hydrogen} functions/atom")
print("\nThe final equation is:")
print(f"{num_carbon} * {functions_per_carbon} + {num_hydrogen} * {functions_per_hydrogen} = {total_functions}")