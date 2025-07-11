# --- Method A ---
# Training an EfficientNet model with 5 insect species
train_A = 36
deploy_A = 13.8
total_A = train_A + deploy_A

# --- Method B ---
# Training an EfficientNet model with 500 insect species
train_B = 126
deploy_B = 13.8
total_B = train_B + deploy_B

# --- Method C ---
# Training a ResNet model with 500 insect species
train_C = 128
deploy_C = 11.8
total_C = train_C + deploy_C

# --- Method D ---
# Manually collecting the data
total_D = 410

print("Calculating the total time for each method:")
print("-" * 40)

print("Method A (5 species):")
print(f"Total time = {train_A} (training) + {deploy_A} (deployment) = {total_A} hours")
print("Result: This method is fast but insufficient as it cannot identify 'all pollinators'.\n")

print("Method B (500 species):")
print(f"Total time = {train_B} (training) + {deploy_B} (deployment) = {total_B} hours")
print("Result: This is a viable and comprehensive method.\n")

print("Method C (500 species):")
print(f"Total time = {train_C} (training) + {deploy_C} (deployment) = {total_C} hours")
print("Result: This is also a viable and comprehensive method.\n")

print("Method D (Manual):")
print(f"Total time = {total_D} (deployment) = {total_D} hours")
print("Result: This method is viable but very time-consuming.\n")

print("-" * 40)
print("Conclusion:")
print(f"Comparing the viable methods (B, C, and D), we find that methods B and C are the fastest at {total_B} hours.")
print("Since both B and C are equally the 'easiest' (i.e., fastest) methods that can achieve the research goal, the correct choice is the one that includes both.")
