import math

# Define the copy number variations (CNVs) and gene weights
cnv_chr1 = 3
cnv_chr2 = 2 # loss of 2
cnv_chr3 = 2

weight_oncogene_A = 0.5
weight_tsg_B = -0.7
weight_oncogene_C = 0.4
# Weights for TSG D, Oncogene E, and TSG F are given for events (loss/gain)
# that do not occur for them, so their score contribution will be zero.

# --- Calculation ---

# Contribution from Chromosome 1 (gain of 3)
# Oncogene A gain contributes positively.
# Gain of Tumor Suppressor D has no specified weight for clonal expansion.
score_chr1 = cnv_chr1 * weight_oncogene_A

# Contribution from Chromosome 2 (loss of 2)
# Loss of Tumor Suppressor B contributes positively. We use the absolute value of its weight.
# Loss of Oncogene E has no specified weight for clonal expansion.
score_chr2 = cnv_chr2 * abs(weight_tsg_B)

# Contribution from Chromosome 3 (gain of 2)
# Oncogene C gain contributes positively.
# Gain of Tumor Suppressor F has no specified weight for clonal expansion.
score_chr3 = cnv_chr3 * weight_oncogene_C

# Total clonal expansion score
total_score = score_chr1 + score_chr2 + score_chr3

# --- Output the results ---
print("Calculating the clonal expansion score based on the provided genetic data:\n")
print(f"Contribution from Oncogene A (Chr 1 Gain): {cnv_chr1} copies * {weight_oncogene_A} weight = {score_chr1}")
print(f"Contribution from Tumor Suppressor B (Chr 2 Loss): {cnv_chr2} copies * {abs(weight_tsg_B)} weight = {score_chr2:.1f}")
print(f"Contribution from Oncogene C (Chr 3 Gain): {cnv_chr3} copies * {weight_oncogene_C} weight = {score_chr3:.1f}\n")
print("Events with no specified weights (gain of a tumor suppressor, loss of an oncogene) contribute 0 to the score.\n")
# Print the final equation with all its components
print(f"Final Equation:")
print(f"{score_chr1} (Oncogene A) + {score_chr2:.1f} (Tumor Suppressor B) + {score_chr3:.1f} (Oncogene C) = {total_score:.1f}\n")
print(f"The final clonal expansion score is: {total_score:.1f}")