# Define the time estimates for each option in hours
time_A_train = 36
time_A_deploy = 13.8

time_B_train = 126
time_B_deploy = 13.8

time_C_train = 128
time_C_deploy = 11.8

time_D_deploy = 410

# Calculate the total time for each option
total_time_A = time_A_train + time_A_deploy
total_time_B = time_B_train + time_B_deploy
total_time_C = time_C_train + time_C_deploy
total_time_D = time_D_deploy

# Print the calculations and results for each option
print("Calculating the total time for each method:")
print(f"Option A (EfficientNet, 5 species): {time_A_train} + {time_A_deploy} = {total_time_A} hours")
print(f"Option B (EfficientNet, 500 species): {time_B_train} + {time_B_deploy} = {total_time_B} hours")
print(f"Option C (ResNet, 500 species): {time_C_train} + {time_C_deploy} = {total_time_C} hours")
print(f"Option D (Manual collection): {total_time_D} hours")

print("\n--- Analysis ---")
print("Option A is insufficient because it only identifies 5 species, not all pollinators.")
print("Comparing the viable options (B, C, and D):")
print(f"Manual collection (D) takes {total_time_D} hours.")
print(f"Both automated methods for 500 species (B and C) take {total_time_B} hours.")
print("Since options B and C take significantly less time than manual collection and are comprehensive enough for the task, they are the easiest methods.")
