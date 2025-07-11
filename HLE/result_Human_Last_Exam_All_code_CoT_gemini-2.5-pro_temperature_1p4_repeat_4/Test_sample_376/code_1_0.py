# Define the copy number changes
gain_chr1 = 3
loss_chr2 = 2
gain_chr3 = 2

# Define the weights for the relevant genes
weight_oncogene_A = 0.5  # per additional copy
weight_suppressor_B = -0.7 # per lost copy
weight_oncogene_C = 0.4  # per additional copy

# Calculate the score contribution from each event
# A gain in an oncogene increases the score.
score_A = gain_chr1 * weight_oncogene_A

# A loss of a tumor suppressor increases the score.
# The effect is positive, so we take the absolute value of the weight or multiply by -1.
score_B = loss_chr2 * (-weight_suppressor_B)

# A gain in an oncogene increases the score.
score_C = gain_chr3 * weight_oncogene_C

# Calculate the total clonal expansion score
total_score = score_A + score_B + score_C

# Print the breakdown of the calculation
print("Clonal Expansion Score Calculation:")
print(f"Oncogene A (Chr 1 Gain): {gain_chr1} copies * {weight_oncogene_A} weight = {score_A}")
print(f"Tumor Suppressor B (Chr 2 Loss): {loss_chr2} copies * {-weight_suppressor_B} effective weight = {score_B}")
print(f"Oncogene C (Chr 3 Gain): {gain_chr3} copies * {weight_oncogene_C} weight = {score_C}")
print(f"Total Score = {score_A} + {score_B} + {score_C} = {total_score}")
