# Time estimates in hours
# Option A: EfficientNet with 5 species
time_train_A = 36
time_deploy_A = 13.8
total_time_A = time_train_A + time_deploy_A

# Option B: EfficientNet with 500 species
time_train_B = 126
time_deploy_B = 13.8
total_time_B = time_train_B + time_deploy_B

# Option C: ResNet with 500 species
time_train_C = 128
time_deploy_C = 11.8
total_time_C = time_train_C + time_deploy_C

# Option D: Manual collection
time_deploy_D = 410
total_time_D = time_deploy_D

print("Analysis of Time Required for Each Method:")
print("-" * 45)

# Print results for Option A
print(f"Option A (EfficientNet, 5 species):")
print(f"Training: {time_train_A} hours + Deployment: {time_deploy_A} hours = Total: {total_time_A} hours")
print("This option is fast but fails to identify all pollinators, so it does not meet the requirements.\n")

# Print results for Option B
print(f"Option B (EfficientNet, 500 species):")
print(f"Training: {time_train_B} hours + Deployment: {time_deploy_B} hours = Total: {total_time_B} hours")
print("This option is comprehensive and meets the requirements.\n")

# Print results for Option C
print(f"Option C (ResNet, 500 species):")
print(f"Training: {time_train_C} hours + Deployment: {time_deploy_C} hours = Total: {total_time_C} hours")
print("This option is also comprehensive and meets the requirements.\n")

# Print results for Option D
print(f"Option D (Manual Collection):")
print(f"Deployment: {time_deploy_D} hours = Total: {total_time_D} hours")
print("This option meets the requirements but is the most time-consuming.\n")

print("-" * 45)
print("Conclusion:")
print(f"Comparing the viable options, methods B and C are the easiest (fastest) with a total time of {total_time_B} hours, compared to {total_time_D} hours for manual collection.")
print("Therefore, both B and C represent the best choices for processing the images efficiently.")
