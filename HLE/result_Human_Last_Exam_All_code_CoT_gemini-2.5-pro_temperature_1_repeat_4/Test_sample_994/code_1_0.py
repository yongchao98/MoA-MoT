# The number of unique four-card combinations with values from 1 to 10.
total_combinations = 715

# The number of those combinations that can be used to make 24.
solvable_combinations = 518

# Calculate the percentage.
percentage = solvable_combinations / total_combinations

# Print the final equation showing how the percentage is calculated.
print("The percentage of solvable combinations is calculated as follows:")
print(f"{solvable_combinations} (solvable) / {total_combinations} (total) = {percentage:.4f}")