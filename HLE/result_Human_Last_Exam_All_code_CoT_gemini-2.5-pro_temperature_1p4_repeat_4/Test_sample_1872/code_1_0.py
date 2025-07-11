# Parameters from the problem statement
# p is the percentage of semantically equivalent but lexically different relations
p = 35

# In a theoretical best-case scenario, the model is limited only by the information it has.
# The JSD (d=0.4) represents a practical difficulty but does not change the theoretical information ceiling.
# The main constraint is the p% of relations that cannot be aligned without external linguistic resources.

# 1. Determine the maximum possible Recall
# The model can, at most, identify the relations that are not lexically different.
# The 35% that are lexically different will be missed, becoming False Negatives.
alignable_fraction = (100 - p) / 100
unalignable_fraction = p / 100
# Recall = TP / (TP + FN) = alignable_fraction / (alignable_fraction + unalignable_fraction)
max_recall = alignable_fraction

# 2. Determine the maximum possible Precision
# A "theoretically maximal" model makes no incorrect alignments, so False Positives = 0.
# Precision = TP / (TP + FP) = alignable_fraction / (alignable_fraction + 0)
max_precision = 1.0

# 3. Calculate the F1 Score using the formula: 2 * (Precision * Recall) / (Precision + Recall)
f1_numerator = 2 * max_precision * max_recall
f1_denominator = max_precision + max_recall
max_f1_score = f1_numerator / f1_denominator

# 4. Print the final result and the equation used to derive it.
print(f"The analysis assumes a perfect model limited only by the p={p}% of lexically different relations, which are fundamentally unalignable without external resources.")
print("")
print(f"Theoretical Maximum Precision = {max_precision}")
print(f"Theoretical Maximum Recall = {max_recall}")
print("")
print("The final F1 score is calculated as follows:")
print(f"F1 Score = (2 * Precision * Recall) / (Precision + Recall)")
print(f"F1 Score = (2 * {max_precision} * {max_recall}) / ({max_precision} + {max_recall}) = {max_f1_score}")