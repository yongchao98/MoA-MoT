# Time estimates for each method in hours
# Option A: Train EfficientNet (5 species)
train_A = 36
deploy_A = 13.8

# Option B: Train EfficientNet (500 species)
train_B = 126
deploy_B = 13.8

# Option C: Train ResNet (500 species)
train_C = 128
deploy_C = 11.8

# Option D: Manual collection
deploy_D = 410

# Calculate total time for each method
total_A = train_A + deploy_A
total_B = train_B + deploy_B
total_C = train_C + deploy_C
total_D = deploy_D # No training time for manual collection

# Print the results
print("Evaluating the total time for each method:")
print(f"Method A (5 species): This method is too limited to identify all pollinators and does not meet the requirements.")
print(f"Total time for Method A = {train_A} + {deploy_A} = {total_A} hours.\n")

print("Method B (500 species):")
print(f"Total time for Method B = {train_B} + {deploy_B} = {total_B} hours.\n")

print("Method C (500 species):")
print(f"Total time for Method C = {train_C} + {deploy_C} = {total_C} hours.\n")

print("Method D (Manual):")
print(f"Total time for Method D = {total_D} hours.\n")

print("Comparison:")
print(f"Manual collection (Method D) takes {total_D} hours.")
print(f"Both machine learning approaches capable of identifying 500 species (Methods B and C) take {total_B} hours.")
print(f"Methods B and C are significantly faster and therefore easier than manual collection (Method D).")
print("Since both B and C are the most efficient and viable methods, the correct choice includes both.")
