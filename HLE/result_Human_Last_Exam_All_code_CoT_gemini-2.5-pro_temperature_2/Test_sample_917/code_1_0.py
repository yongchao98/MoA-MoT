# Option A: EfficientNet with 5 species
train_A = 36
deploy_A = 13.8
total_A = train_A + deploy_A

# Option B: EfficientNet with 500 species
train_B = 126
deploy_B = 13.8
total_B = train_B + deploy_B

# Option C: ResNet with 500 species
train_C = 128
deploy_C = 11.8
total_C = train_C + deploy_C

# Option D: Manual collection
total_D = 410

print("Time Analysis for Pollinator Identification Methods:")
print("-" * 50)
print(f"Option A Total Time: {train_A} + {deploy_A} = {total_A} hours (Note: Fails to meet goal of identifying all pollinators)")
print(f"Option B Total Time: {train_B} + {deploy_B} = {total_B} hours")
print(f"Option C Total Time: {train_C} + {deploy_C} = {total_C} hours")
print(f"Option D Total Time: {total_D} hours")
print("-" * 50)
print("Conclusion:")
print("Methods B and C are the most effective and time-efficient options. They both take the same total time and are significantly faster than manual collection.")
