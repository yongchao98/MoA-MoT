# Toluene molecule formula: C7H8
num_C = 7
num_H = 8

# Number of contracted Gaussian functions for C and H in the 6-31G basis set
# For Carbon (a heavy atom):
# 1s (core) -> 1 contracted function
# 2s, 2px, 2py, 2pz (inner valence) -> 4 contracted functions
# 2s, 2px, 2py, 2pz (outer valence) -> 4 contracted functions
# Total per C = 1 + 4 + 4 = 9
cgf_per_C = 9

# For Hydrogen:
# 1s (inner valence) -> 1 contracted function
# 1s (outer valence) -> 1 contracted function
# Total per H = 1 + 1 = 2
cgf_per_H = 2

# Calculate the total number of contracted functions for each element
total_C = num_C * cgf_per_C
total_H = num_H * cgf_per_H

# Calculate the grand total
total_cgf = total_C + total_H

# Print the final result and the equation
print("Calculation for Toluene (C7H8) with 6-31G basis set:")
print(f"Total from Carbon: {num_C} atoms * {cgf_per_C} functions/atom = {total_C}")
print(f"Total from Hydrogen: {num_H} atoms * {cgf_per_H} functions/atom = {total_H}")
print(f"Total contracted Gaussian functions = {total_C} + {total_H} = {total_cgf}")