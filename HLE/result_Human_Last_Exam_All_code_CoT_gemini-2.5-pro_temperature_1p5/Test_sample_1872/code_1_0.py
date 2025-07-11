# Define the percentage of semantically equivalent but lexically different relations.
p = 35

# Convert the percentage to a decimal.
p_decimal = p / 100.0

# In the absence of external linguistic resources, the p% of lexically different relations
# are fundamentally unalignable. This creates a hard limit on recall.
# 1. We calculate the maximum theoretical recall. In the best case, the model identifies
# all relations except for the unalignable ones.
max_recall = 1.0 - p_decimal

# 2. To achieve the theoretical maximum F1 score, we assume the model makes no mistakes
# on the alignable portion, meaning it has no false positives. This results in a
# precision of 1.0.
max_precision = 1.0

# 3. We calculate the maximal F1 score using the formula:
# F1 = 2 * (Precision * Recall) / (Precision + Recall)
max_f1_score = (2 * (max_precision * max_recall)) / (max_precision + max_recall)

# Print the final result along with the values used in the equation.
print("This problem assesses the theoretical limits of knowledge graph alignment without external linguistic aid.")
print(f"Given that {p}% of relations are lexically different, they are considered unalignable.")
print("\nStep 1: Calculate Maximum Theoretical Precision")
print(f"In a perfect scenario, no incorrect alignments are made from the alignable portion.")
print(f"Maximal Precision = {max_precision}")

print("\nStep 2: Calculate Maximum Theoretical Recall")
print(f"The {p}% of unalignable relations cannot be found.")
print(f"Maximal Recall = 1.0 - {p_decimal} = {max_recall}")

print("\nStep 3: Calculate Maximal F1 Score")
print(f"F1 Score = 2 * (Precision * Recall) / (Precision + Recall)")
print(f"Final Equation: {max_f1_score:.5f} = 2 * ({max_precision} * {max_recall}) / ({max_precision} + {max_recall})")

print(f"\nThe theoretically maximal F1 score is: {max_f1_score:.5f}")
<<<0.78788>>>