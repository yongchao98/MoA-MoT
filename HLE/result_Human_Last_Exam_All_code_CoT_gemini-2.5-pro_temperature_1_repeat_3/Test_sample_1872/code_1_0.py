# p: The percentage of relations that are semantically equivalent but lexically different.
p = 35

# Convert the percentage to a proportion.
p_prop = p / 100.0

# In a "theoretically maximal" scenario, the model makes no wrong guesses,
# so False Positives (FP) = 0. This makes Precision = 1.
# Precision = TP / (TP + FP)
precision = 1.0

# The model can only identify the lexically identical relations. The proportion of these
# is (1 - p_prop). These are the True Positives (TP).
# The model misses the lexically different but semantically equivalent relations.
# The proportion of these is p_prop. These are the False Negatives (FN).
# Recall = TP / (TP + FN) = (1 - p_prop) / ((1 - p_prop) + p_prop) = 1 - p_prop
recall = 1.0 - p_prop

# The F1 score is the harmonic mean of Precision and Recall.
# F1 = 2 * (Precision * Recall) / (Precision + Recall)
f1_score = 2 * (precision * recall) / (precision + recall)

# Print the final calculation as requested, showing each number in the equation.
print(f"Theoretically Maximal F1 Score Calculation:")
print(f"p = {p}%")
print(f"Precision = {precision:.2f}")
print(f"Recall = 1.0 - ({p}/100.0) = {recall:.2f}")
print(f"F1 Score = 2 * ({precision:.2f} * {recall:.2f}) / ({precision:.2f} + {recall:.2f}) = {f1_score}")