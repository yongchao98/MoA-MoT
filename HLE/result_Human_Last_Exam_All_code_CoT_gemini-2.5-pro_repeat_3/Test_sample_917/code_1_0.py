# Time estimates for each method in hours
# Method A: EfficientNet with 5 species
time_A_train = 36
time_A_deploy = 13.8

# Method B: EfficientNet with 500 species
time_B_train = 126
time_B_deploy = 13.8

# Method C: ResNet with 500 species
time_C_train = 128
time_C_deploy = 11.8

# Method D: Manual data collection
time_D_deploy = 410

# Calculate total time for each method
total_time_A = time_A_train + time_A_deploy
total_time_B = time_B_train + time_B_deploy
total_time_C = time_C_train + time_C_deploy
total_time_D = time_D_deploy # No training time for manual collection

# Print the breakdown and total time for each method
print("Calculating the total time for each method:")

# Method A is insufficient for the task, but we calculate its time for completeness.
print(f"Method A: {time_A_train} (training) + {time_A_deploy} (deployment) = {total_time_A} hours")

# Method B is a viable automated approach.
print(f"Method B: {time_B_train} (training) + {time_B_deploy} (deployment) = {total_time_B} hours")

# Method C is another viable automated approach.
print(f"Method C: {time_C_train} (training) + {time_C_deploy} (deployment) = {total_time_C} hours")

# Method D is the manual approach.
print(f"Method D: 0 (training) + {time_D_deploy} (deployment) = {total_time_D} hours")

print("\nComparison:")
print("Method A is the fastest but does not meet the scientific goal of identifying all pollinators.")
print("Comparing the viable methods (B, C, and D):")
print(f"Method B takes {total_time_B} hours.")
print(f"Method C takes {total_time_C} hours.")
print(f"Method D takes {total_time_D} hours.")
print("Methods B and C require the least amount of time and are equally easy (fast). Both are significantly easier than manual collection.")
