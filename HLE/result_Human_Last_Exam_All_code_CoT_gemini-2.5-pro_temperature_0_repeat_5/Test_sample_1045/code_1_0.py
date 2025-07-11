# Number of atoms in toluene (C7H8)
num_C = 7
num_H = 8

# Number of contracted Gaussian functions per atom for the 6-31G basis set
funcs_per_C = 9  # 1 (core) + 4 (inner valence) + 4 (outer valence)
funcs_per_H = 2  # 1 (inner valence) + 1 (outer valence)

# Calculate the total number of functions from Carbon atoms
total_C = num_C * funcs_per_C

# Calculate the total number of functions from Hydrogen atoms
total_H = num_H * funcs_per_H

# Calculate the grand total for the molecule
grand_total = total_C + total_H

# Print the explanation and the final calculation
print("For a 6-31G basis set:")
print(f"- Each Carbon atom (C) has {funcs_per_C} contracted basis functions.")
print(f"- Each Hydrogen atom (H) has {funcs_per_H} contracted basis functions.")
print("\nFor Toluene (C7H8):")
print(f"Total functions = ({num_C} C atoms * {funcs_per_C} functions/C) + ({num_H} H atoms * {funcs_per_H} functions/H)")
print(f"                = {total_C} + {total_H}")
print(f"                = {grand_total}")

print(f"\nTherefore, a total of {grand_total} contracted Gaussian functions are used.")