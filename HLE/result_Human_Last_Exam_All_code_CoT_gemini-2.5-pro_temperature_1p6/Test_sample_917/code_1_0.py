# Step 1: Define the time estimates for each method.
# Times are in hours. Format is (preparation_time, deployment_time).
times = {
    'A': (36, 13.8),
    'B': (126, 13.8),
    'C': (128, 11.8),
    'D': (0, 410) # Manual collection has no model training time.
}

# Step 2 & 3: Calculate and print the total time for each method.
print("--- Calculating Total Time for Each Method ---\n")

# Method A
prep_A, deploy_A = times['A']
total_A = prep_A + deploy_A
print(f"Method A (EfficientNet, 5 species):")
print(f"Equation: {prep_A} + {deploy_A} = {total_A} hours\n")

# Method B
prep_B, deploy_B = times['B']
total_B = prep_B + deploy_B
print(f"Method B (EfficientNet, 500 species):")
print(f"Equation: {prep_B} + {deploy_B} = {total_B} hours\n")

# Method C
prep_C, deploy_C = times['C']
total_C = prep_C + deploy_C
print(f"Method C (ResNet, 500 species):")
print(f"Equation: {prep_C} + {deploy_C} = {total_C} hours\n")

# Method D
prep_D, deploy_D = times['D']
total_D = prep_D + deploy_D
print(f"Method D (Manual Collection):")
print(f"Equation: {prep_D} + {deploy_D} = {total_D} hours\n")


# Step 4: Analyze the results to find the easiest method.
print("--- Analysis ---\n")
print("To find the easiest method, we need the fastest option that meets the requirement of identifying all pollinators.")
print("1. Method A is very fast ({:.1f} hours), but training on only 5 species is insufficient for identifying 'all pollinators' in a diverse ecosystem.".format(total_A))
print("2. Method D ({:.1f} hours) is thorough but is the most time-consuming and therefore the hardest method.".format(total_D))
print("3. Methods B ({:.1f} hours) and C ({:.1f} hours) are both trained on a comprehensive 500 species.".format(total_B, total_C))
print("\nComparing the viable options, Methods B and C take the least amount of time ({:.1f} hours) to complete the task properly.".format(total_B))
print("Since both B and C represent the easiest and most effective approach, the correct answer choice encompasses both.")
print("\nFinal Conclusion: The easiest methods are B and C.")
