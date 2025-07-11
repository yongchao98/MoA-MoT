# p is the percentage of semantically equivalent but lexically different relations.
p = 35

# In the theoretical best-case scenario without external linguistic resources,
# a model can only align relations that are not lexically different.
# It is assumed that the model makes no mistakes on this alignable subset.

# Precision = TP / (TP + FP)
# A perfect model makes no incorrect alignments, so False Positives (FP) = 0.
# This makes Precision = TP / TP = 1.0.
precision = 1.0

# Recall = TP / (TP + FN)
# The model can only identify the relations that are NOT lexically different.
# The percentage of these relations is (100 - p)%.
# So, the maximum theoretical recall is (100 - p) / 100.
recall = (100 - p) / 100

# F1 Score is the harmonic mean of Precision and Recall.
# F1 = 2 * (Precision * Recall) / (Precision + Recall)
f1_score = 2 * (precision * recall) / (precision + recall)

# --- Output the calculation steps ---
print(f"Given p = {p}%\n")

print("Step 1: Determine the theoretical maximum Precision.")
print("A perfect model makes no incorrect predictions, so FP = 0, and Precision = 1.0\n")

print("Step 2: Determine the theoretical maximum Recall.")
print(f"The model can only align the (100 - p)% of relations that are not lexically different.")
print(f"Recall = (100 - {p}) / 100 = {recall}\n")

print("Step 3: Calculate the F1 score.")
print("F1 Score = 2 * (Precision * Recall) / (Precision + Recall)")
print(f"F1 Score = 2 * ({precision} * {recall}) / ({precision} + {recall})")
print(f"F1 Score = {2 * precision * recall} / {precision + recall}")
print(f"Final F1 Score = {f1_score}")