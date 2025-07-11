# Define the copy number changes from the problem description
gain_chr1 = 3
loss_chr2 = 2
gain_chr3 = 2

# Define the weights for the relevant genes
weight_oncogene_A = 0.5  # per additional copy
weight_tsg_B_per_loss = -0.7  # per lost copy
weight_oncogene_C = 0.4  # per additional copy

# --- Calculation ---

# 1. Oncogene A is on Chromosome 1 (gain of 3 copies)
# Contribution = number of copies gained * weight per copy
score_A = gain_chr1 * weight_oncogene_A

# 2. Tumor Suppressor B is on Chromosome 2 (loss of 2 copies)
# Losing a tumor suppressor promotes expansion, so its contribution to the score is positive.
# We take the negative of the provided weight to find its positive impact on the score.
score_B = loss_chr2 * (-weight_tsg_B_per_loss)

# 3. Oncogene C is on Chromosome 3 (gain of 2 copies)
# Contribution = number of copies gained * weight per copy
score_C = gain_chr3 * weight_oncogene_C

# The other genes (D, E, F) do not contribute, as the type of change
# (gain of a tumor suppressor or loss of an oncogene) does not have a defined weight.

# 4. Calculate the total score
total_score = score_A + score_B + score_C

# --- Output ---

print("Calculating the clonal expansion score based on the relevant genetic events:")
print("-" * 60)

# Explanation for Oncogene A
print(f"Oncogene A (Chromosome 1 Gain): {gain_chr1} copies * {weight_oncogene_A} score/copy = {score_A:.1f}")

# Explanation for Tumor Suppressor B
print(f"Tumor Suppressor B (Chromosome 2 Loss): {loss_chr2} copies * {-weight_tsg_B_per_loss:.1f} score/copy = {score_B:.1f}")

# Explanation for Oncogene C
print(f"Oncogene C (Chromosome 3 Gain): {gain_chr3} copies * {weight_oncogene_C} score/copy = {score_C:.1f}")
print("-" * 60)

# Print the final equation with all numbers
print("Total Clonal Expansion Score:")
print(f"{score_A:.1f} (from Oncogene A) + {score_B:.1f} (from Tumor Suppressor B) + {score_C:.1f} (from Oncogene C) = {total_score:.1f}")
print("-" * 60)
