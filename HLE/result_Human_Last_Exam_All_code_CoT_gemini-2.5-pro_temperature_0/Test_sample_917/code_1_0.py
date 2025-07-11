# Time estimates in hours
time_A_train = 36
time_A_deploy = 13.8

time_B_train = 126
time_B_deploy = 13.8

time_C_train = 128
time_C_deploy = 11.8

time_D_deploy = 410

# Calculate total time for each option
total_time_A = time_A_train + time_A_deploy
total_time_B = time_B_train + time_B_deploy
total_time_C = time_C_train + time_C_deploy
total_time_D = time_D_deploy

# Print the results
print(f"Option A (EfficientNet, 5 species) Total Time: {time_A_train} + {time_A_deploy} = {total_time_A} hours")
print("Note: This option is inadequate as it fails to identify the majority of pollinator species.")
print("-" * 20)
print(f"Option B (EfficientNet, 500 species) Total Time: {time_B_train} + {time_B_deploy} = {total_time_B} hours")
print(f"Option C (ResNet, 500 species) Total Time: {time_C_train} + {time_C_deploy} = {total_time_C} hours")
print(f"Option D (Manual) Total Time: {time_D_deploy} = {total_time_D} hours")
print("-" * 20)
print("Conclusion: Options B and C are the most time-efficient methods that meet the project's requirements. They both take the same amount of total time.")
