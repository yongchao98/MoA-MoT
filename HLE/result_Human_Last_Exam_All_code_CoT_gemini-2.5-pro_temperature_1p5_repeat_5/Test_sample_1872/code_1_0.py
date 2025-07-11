import sys

# Step 1: Define the given parameters
p_percent = 35  # percentage of semantically equivalent but lexically different relations
d = 0.4         # Jensen-Shannon divergence between relational distributions

# Convert percentage to a proportion
p = p_percent / 100.0

# Step 2: Set Precision for the theoretical maximum
# To achieve the maximal F1 score, we assume the model makes no incorrect alignments.
# This means False Positives (FP) = 0, which results in a Precision of 1.0.
precision = 1.0
print(f"To calculate the theoretical maximum F1 score, we assume perfect precision.")
print(f"Precision = 1.0\n")

# Step 3: Calculate the maximum possible Recall
print(f"Next, we calculate the maximum achievable Recall.")
# Total relations are partitioned into two groups.
# Group 1: (1 - p) are easily alignable.
# Group 2: p are lexically different and can only be aligned based on structure.
recall_from_easy_group = (1 - p)
print(f"The proportion of easily alignable relations is (100% - p) = {100 - p_percent}%. Their contribution to recall is {recall_from_easy_group:.2f}.")

# For the lexically different group, the success rate is (1 - d).
recall_from_hard_group = p * (1 - d)
print(f"The proportion of lexically different relations is p = {p_percent}%.")
print(f"The structural dissimilarity is d = {d:.1f}, so the fraction that can be aligned is (1 - d) = {1 - d:.1f}.")
print(f"Their contribution to recall is p * (1 - d) = {p:.2f} * {1 - d:.1f} = {recall_from_hard_group:.2f}.")

# Total Recall is the sum of recall from both groups.
recall = recall_from_easy_group + recall_from_hard_group
print(f"\nTotal Recall = (Recall from easy group) + (Recall from hard group)")
print(f"Total Recall = {recall_from_easy_group:.2f} + {recall_from_hard_group:.2f} = {recall:.2f}\n")


# Step 4: Calculate the F1 Score
# F1 = 2 * (Precision * Recall) / (Precision + Recall)
if (precision + recall) == 0:
    f1_score = 0.0
else:
    f1_score = 2 * (precision * recall) / (precision + recall)

# Output the final calculation clearly
print("Finally, we calculate the F1 score using the formula: F1 = 2 * (Precision * Recall) / (Precision + Recall)")
print(f"F1 = 2 * ({precision:.1f} * {recall:.2f}) / ({precision:.1f} + {recall:.2f})")
numerator = 2 * precision * recall
denominator = precision + recall
print(f"F1 = {numerator:.2f} / {denominator:.2f}")
print(f"F1 Score = {f1_score:.4f}")

# Suppress all previous output and print only the final answer for the platform.
# We redirect stdout, print the final answer, and then restore stdout.
original_stdout = sys.stdout
sys.stdout = open('/dev/null', 'w')
print(f"\n<<<{(2 * (1.0 * ( (1-p) + p*(1-d) )) / (1.0 + ( (1-p) + p*(1-d) ))):.4f}>>>", file=original_stdout, end='')
sys.stdout.close()
sys.stdout = original_stdout