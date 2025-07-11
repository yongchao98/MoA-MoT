# The percentage of semantically equivalent but lexically different relations.
p = 35

# In a theoretical best-case scenario, the model is assumed to be perfectly precise,
# meaning it makes no incorrect alignments (False Positives = 0).
# Precision = True Positives / (True Positives + False Positives) = 1
precision = 1.0

# The model's recall is limited by the information it has. Without external
# linguistic resources, it cannot align the p% of lexically different relations.
# These become unavoidable False Negatives. The maximum possible recall is the
# fraction of relations that are alignable.
# Recall = True Positives / (True Positives + False Negatives) = (100 - p) / 100
recall = (100.0 - p) / 100.0

# The F1 score is the harmonic mean of Precision and Recall.
f1_score = 2 * (precision * recall) / (precision + recall)

# Print the final equation with the calculated numbers to show the steps.
print("To find the theoretically maximal F1 score, we make the following assumptions:")
print("1. A perfect model makes no incorrect alignments, so Precision = 1.")
print(f"2. With p={p}%, the {p}% of lexically different relations cannot be aligned, limiting the maximal Recall.")
print("\nCalculation:")
print("F1 = 2 * (Precision * Recall) / (Precision + Recall)")
print(f"F1 = 2 * ({precision} * {recall}) / ({precision} + {recall})")
print(f"F1 = {2 * precision * recall} / {precision + recall}")
print(f"F1 = {f1_score}")