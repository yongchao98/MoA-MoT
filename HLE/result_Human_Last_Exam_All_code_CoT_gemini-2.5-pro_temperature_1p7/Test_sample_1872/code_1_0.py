import sys

# The percentage of semantically equivalent but lexically different relations.
p = 35

# Convert the percentage to a proportion.
# This represents the proportion of true alignments that are impossible to find.
# These will be the False Negatives (FN) in a best-case scenario.
unalignable_proportion = p / 100.0

# In a theoretically maximal scenario, we assume an ideal model.
# 1. It makes no incorrect alignments. So, False Positives (FP) = 0.
#    This leads to a Precision of 1.0.
#    Precision = TP / (TP + FP) = TP / (TP + 0) = 1.0
precision = 1.0

# 2. It correctly aligns all relations that are not lexically different.
#    The proportion of alignable relations (True Positives, TP) is 1 - unalignable_proportion.
#    Recall = TP / (TP + FN).
#    This simplifies to the proportion of relations that the model can successfully align.
recall = 1.0 - unalignable_proportion

# F1 score is the harmonic mean of precision and recall.
# F1 = 2 * (Precision * Recall) / (Precision + Recall)
f1_score = (2 * precision * recall) / (precision + recall)

# To fulfill the request of printing each number in the final equation.
numerator = 2 * precision * recall
denominator = precision + recall

# Print the explanation and the step-by-step calculation.
print("This script calculates the theoretically maximal F1 score under the given constraints.")
print(f"\nGiven p = {p}%, which is the percentage of relations that are lexically different.")
print("A model without external linguistic resources cannot align these relations.")
print("\nWe assume an ideal model where:")
print(f" - False Positives = 0, leading to Precision = {precision}")
print(f" - All alignable relations are found, leading to Recall = 1 - ({p}/100) = {recall:.2f}")

print("\nThe F1 Score is calculated as follows:")
print(f"F1 Score = 2 * (Precision * Recall) / (Precision + Recall)")
print(f"F1 Score = 2 * ({precision} * {recall:.2f}) / ({precision} + {recall:.2f})")
print(f"F1 Score = {numerator:.2f} / {denominator:.2f}")
print(f"Maximal F1 Score = {f1_score:.4f}")