import sys

# p: The percentage of semantically equivalent but lexically different relations.
p = 35

# Convert the percentage to a decimal fraction.
p_frac = p / 100

# For a theoretically perfect model under these constraints:
# True Positives (TP) are the (1 - p_frac) portion of relations that are discoverable.
# False Positives (FP) are 0, as a perfect model makes no mistakes.
# False Negatives (FN) are the p_frac portion of relations that are undiscoverable.

# Precision = TP / (TP + FP). Since FP=0, Precision is 1.
precision = 1.0

# Recall = TP / (TP + FN). This is the fraction of all true alignments that are discoverable.
recall = 1.0 - p_frac

# The F1 score is the harmonic mean of Precision and Recall.
# F1 = 2 * (Precision * Recall) / (Precision + Recall)
f1_score = 2 * (precision * recall) / (precision + recall)

# --- Output ---
# The problem asks to output the final equation with each number.

# First, print the derived values for precision and recall.
print(f"Given p = {p}% of lexically different relations:")
print(f"The theoretical maximum Precision is {precision}.")
print(f"The theoretical maximum Recall is 1.0 - {p_frac} = {recall:.2f}.")
print("-" * 30)

# Now, print the F1 score calculation step-by-step with the numbers.
print("The maximal F1 score is calculated as follows:")
print(f"F1 = 2 * (Precision * Recall) / (Precision + Recall)")
# Print the equation with the numeric values substituted.
print(f"F1 = 2 * ({precision} * {recall:.2f}) / ({precision} + {recall:.2f})")
# Print the evaluation of the numerator and denominator.
numerator = 2 * precision * recall
denominator = precision + recall
print(f"F1 = {numerator:.2f} / {denominator:.2f}")
# Print the final result.
print(f"F1 = {f1_score}")

# The final answer format for the platform.
# Redirecting to stdout to be captured by the platform.
sys.stdout = sys.__stderr__
print(f'<<<{f1_score}>>>')
