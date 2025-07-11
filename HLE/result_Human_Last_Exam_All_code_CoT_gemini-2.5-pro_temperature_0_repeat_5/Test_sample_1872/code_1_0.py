# Given parameters
# p: percentage of semantically equivalent but lexically different relations
p = 35
# d: Jensen-Shannon divergence between the relational distributions
d = 0.4

# Convert percentage to a fraction
p_frac = p / 100.0

print(f"Step 1: Define the model for maximum recall.")
print(f"The relations are split into two groups:")
print(f"- {100-p}% are lexically similar, assumed to be perfectly alignable (Recall = 1.0).")
print(f"- {p}% are lexically different and rely on structure for alignment.")
print(f"The success rate for structure-based alignment is limited by the divergence d, modeled as (1 - d).")
print(f"Max Recall = (1 - p/100) * 1.0 + (p/100) * (1 - d)\n")

# Calculate the maximum possible recall
# Contribution from lexically similar relations
recall_lex_similar = (1 - p_frac) * 1.0
# Contribution from lexically different relations
recall_lex_different = p_frac * (1 - d)
# Total maximum recall
max_recall = recall_lex_similar + recall_lex_different

print(f"Step 2: Calculate the maximum recall with p = {p} and d = {d}.")
print(f"Max Recall = (1 - {p_frac}) + {p_frac} * (1 - {d})")
print(f"Max Recall = {1-p_frac} + {p_frac} * {1-d}")
print(f"Max Recall = {1-p_frac} + {p_frac * (1-d)}")
print(f"Max Recall = {max_recall}\n")

# To achieve the maximal F1 score, we assume the model makes no mistakes,
# meaning Precision = 1.0.
precision = 1.0

# Calculate the maximum F1 score
# F1 = 2 * (Precision * Recall) / (Precision + Recall)
max_f1 = 2 * (precision * max_recall) / (precision + max_recall)

print(f"Step 3: Calculate the maximum F1 score.")
print(f"Assuming Precision = {precision} to maximize the F1 score.")
print(f"Max F1 = 2 * (Precision * Recall) / (Precision + Recall)")
print(f"Max F1 = 2 * ({precision} * {max_recall}) / ({precision} + {max_recall})")
print(f"Max F1 = {2 * precision * max_recall} / {precision + max_recall}")
print(f"Max F1 = {max_f1}")

# Final answer
# print(f"\n<<<The theoretically maximal F1 score is {max_f1:.4f}>>>")
<<<0.9247>>>