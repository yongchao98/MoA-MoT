# Time estimates in hours for each method
# A: EfficientNet with 5 species
time_train_A = 36
time_deploy_A = 13.8

# B: EfficientNet with 500 species
time_train_B = 126
time_deploy_B = 13.8

# C: ResNet with 500 species
time_train_C = 128
time_deploy_C = 11.8

# D: Manual data collection
time_deploy_D = 410

# --- Calculations ---
# Total time for option A
total_time_A = time_train_A + time_deploy_A
# Total time for option B
total_time_B = time_train_B + time_deploy_B
# Total time for option C
total_time_C = time_train_C + time_deploy_C
# Total time for option D
total_time_D = time_deploy_D

print("Analysis of Time Commitment for Each Method:")
print("-" * 45)
# Print results for A, noting its limitation
print(f"A. EfficientNet (5 species): {time_train_A} + {time_deploy_A} = {total_time_A} hours")
print("   (Note: This is insufficient as it only identifies 5 species.)\n")

# Print results for B
print(f"B. EfficientNet (500 species): {time_train_B} + {time_deploy_B} = {total_time_B} hours\n")

# Print results for C
print(f"C. ResNet (500 species): {time_train_C} + {time_deploy_C} = {total_time_C} hours\n")

# Print results for D
print(f"D. Manual Collection: {total_time_D} hours\n")

print("-" * 45)
print("Conclusion:")
print("Methods B and C are the easiest as they take the least amount of time (139.8 hours) while being capable of identifying a wide range of species.")
print("Since both are equally efficient, the best answer choice is F, which includes both B and C.")
