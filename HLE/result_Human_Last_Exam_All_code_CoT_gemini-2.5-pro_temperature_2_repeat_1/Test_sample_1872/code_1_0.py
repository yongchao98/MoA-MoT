import math

# --- Given Parameters ---
# p: percentage of semantically equivalent but lexically different relations
p_percent = 35
p = p_percent / 100.0

# d: Jensen-Shannon divergence representing structural differences
d = 0.4

# --- Calculation ---

# Step 1: Calculate the proportion of relations that are lexically alignable.
# Since p% are lexically different and we have no external resources, they cannot be aligned.
# The remaining portion is potentially alignable.
lexically_alignable_proportion = 1 - p

# Step 2: Of the lexically alignable relations, calculate the proportion that are structurally alignable.
# The divergence 'd' makes this fraction of the relations structurally unalignable.
structurally_alignable_proportion = 1 - d

# Step 3: Calculate the maximum theoretical Recall.
# Recall is the fraction of total true alignments that we can correctly identify.
# This is the product of the proportions that survive both the lexical and structural filters.
max_recall = lexically_alignable_proportion * structurally_alignable_proportion

# Step 4: Determine the theoretical maximum Precision.
# For a "theoretically maximal" F1 score, we assume an ideal model that achieves perfect
# precision by never making an incorrect alignment (False Positives = 0).
max_precision = 1.0

# Step 5: Calculate the maximal F1 score.
# F1 = 2 * (Precision * Recall) / (Precision + Recall)
max_f1_score = 2 * (max_precision * max_recall) / (max_precision + max_recall)


# --- Output the results ---

print(f"Given Parameters:")
print(f"  - Percentage of lexically different relations (p): {p_percent}%")
print(f"  - Jensen-Shannon divergence (d): {d}\n")

print("--- Calculation Steps ---")

print(f"1. Proportion of lexically alignable relations = 1 - p")
print(f"   = 1 - {p} = {lexically_alignable_proportion}\n")

print(f"2. Proportion of structurally alignable relations = 1 - d")
print(f"   = 1 - {d} = {structurally_alignable_proportion}\n")

print(f"3. Maximum theoretical Recall = (1 - p) * (1 - d)")
print(f"   = ({lexically_alignable_proportion}) * ({structurally_alignable_proportion}) = {max_recall:.4f}\n")

print(f"4. Maximum theoretical Precision is assumed to be 1.0 for a maximal F1 score.\n")

print(f"5. Maximum theoretical F1 Score = 2 * (Precision * Recall) / (Precision + Recall)")
print(f"   = 2 * ({max_precision} * {max_recall:.4f}) / ({max_precision} + {max_recall:.4f})")
print(f"   = {2 * max_precision * max_recall:.4f} / {max_precision + max_recall:.4f}")
print(f"   = {max_f1_score:.4f}\n")

print(f"The theoretically maximal F1 score is {max_f1_score}.")