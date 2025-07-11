# Parameters from the problem description
# p: percentage of semantically equivalent but lexically different relations
p_percent = 35
# d: Jensen-Shannon divergence between relational distributions
d = 0.4

# Convert percentage to a proportion
p = p_percent / 100.0

# In a "theoretically maximal" F1 score scenario, we assume the model is tuned
# to have perfect precision (P=1.0) to avoid any false positives.
# The F1 score is then determined solely by the maximal achievable recall (R).
# F1 = 2 * P * R / (P + R). With P=1, F1 = 2 * R / (1 + R).
precision = 1.0

# The set of true alignments can be split into two groups:
# 1. Lexically similar relations: The model can identify these perfectly.
#    Their proportion is (1 - p).
prop_lexically_similar = 1 - p
recall_lexically_similar = 1.0

# 2. Lexically different relations: These must be aligned using graph structure alone.
#    The success rate of structural alignment is modeled as (1 - d).
#    Their proportion is p.
prop_lexically_different = p
recall_structural = 1 - d

# The total recall is the weighted sum of the recall from each group.
# R = (proportion_of_group_1 * recall_of_group_1) + (proportion_of_group_2 * recall_of_group_2)
recall = (prop_lexically_similar * recall_lexically_similar) + (prop_lexically_different * recall_structural)

# Calculate the final F1 score
f1_score = (2 * precision * recall) / (precision + recall)

# --- Output the results and the formula steps ---
print("Step 1: Calculate the total Recall (R)")
print(f"The overall recall is a weighted average based on the type of relation.")
print(f"   - For the {(1-p):.2f} portion of lexically similar relations, assume Recall = {recall_lexically_similar}")
print(f"   - For the {p:.2f} portion of lexically different relations, success depends on structure (Recall = 1 - d = {recall_structural:.2f})")
print("\nEquation for Recall (R):")
# Here we print the final equation with all the numbers.
print(f"R = ({1-p:.2f} * {recall_lexically_similar:.1f}) + ({p:.2f} * (1 - {d:.1f}))")
print(f"R = {prop_lexically_similar:.2f} + {prop_lexically_different:.2f} * {recall_structural:.2f}")
print(f"R = {prop_lexically_similar:.2f} + {prop_lexically_different * recall_structural:.2f}")
print(f"R = {recall:.2f}")

print("\nStep 2: Calculate the F1 Score")
print("Assuming theoretical maximum where Precision (P) = 1.0")
print("\nEquation for F1 Score:")
# Here we print the final equation with all the numbers.
print(f"F1 = (2 * P * R) / (P + R)")
print(f"F1 = (2 * {precision:.1f} * {recall:.2f}) / ({precision:.1f} + {recall:.2f})")
print(f"F1 = {2 * precision * recall:.2f} / {precision + recall:.2f}")
print(f"F1 = {f1_score}")

# Final answer to be parsed
# print(f"<<<{f1_score:.4f}>>>")