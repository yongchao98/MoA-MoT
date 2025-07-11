# Define the copy number changes and weights based on the problem description.

# --- Chromosome 1 ---
# Gain of 3 copies, affecting Oncogene A.
copies_gained_A = 3
weight_A = 0.5 # per additional copy

# --- Chromosome 2 ---
# Loss of 2 copies, affecting Tumor Suppressor B.
copies_lost_B = 2
# The weight is -0.7 per lost copy. A loss increases clonal expansion, so its effect is positive.
# The score contribution is copies_lost * (-weight), which is 2 * -(-0.7) = 2 * 0.7.
effective_weight_B = 0.7

# --- Chromosome 3 ---
# Gain of 2 copies, affecting Oncogene C.
copies_gained_C = 2
weight_C = 0.4 # per additional copy

# Calculate the score contribution from each relevant gene
score_A = copies_gained_A * weight_A
score_B = copies_lost_B * effective_weight_B
score_C = copies_gained_C * weight_C

# Calculate the total clonal expansion score
total_score = score_A + score_B + score_C

# Print the final equation and the result
# The problem asks to show each number in the final equation.
print("Clonal Expansion Score Calculation:")
print(f"Contribution from Oncogene A (gain): {copies_gained_A} copies * {weight_A} weight = {score_A}")
print(f"Contribution from Tumor Suppressor B (loss): {copies_lost_B} copies * {effective_weight_B} weight = {score_B:.1f}")
print(f"Contribution from Oncogene C (gain): {copies_gained_C} copies * {weight_C} weight = {score_C}")
print("-" * 35)
print(f"Total Score = ({copies_gained_A} * {weight_A}) + ({copies_lost_B} * {effective_weight_B}) + ({copies_gained_C} * {weight_C}) = {total_score}")