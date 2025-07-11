# Define the weights and copy number changes for each relevant gene
oncogene_A_gain = 3
oncogene_A_weight = 0.5

tsg_B_loss = 2
tsg_B_weight = 0.7  # Loss of a tumor suppressor promotes expansion, so we use a positive value

oncogene_C_gain = 2
oncogene_C_weight = 0.4

# The repressor causes a functional loss of Tumor Suppressor D.
# We assume this is equivalent to the loss of the 2 baseline copies.
tsg_D_functional_loss = 2
tsg_D_weight = 0.6 # Loss of a tumor suppressor promotes expansion

# Calculate the contribution of each event to the clonal expansion score
score_A = oncogene_A_gain * oncogene_A_weight
score_B = tsg_B_loss * tsg_B_weight
score_C = oncogene_C_gain * oncogene_C_weight
score_D = tsg_D_functional_loss * tsg_D_weight

# Calculate the total clonal expansion score
total_score = score_A + score_B + score_C + score_D

# Print the breakdown of the calculation and the final score
# The final equation shows each number that contributes to the total score.
print("Clonal Expansion Score Calculation:")
print(f"Oncogene A (Gain): {oncogene_A_gain} copies * {oncogene_A_weight}/copy = {score_A}")
print(f"Tumor Suppressor B (Loss): {tsg_B_loss} copies * {tsg_B_weight}/copy = {score_B}")
print(f"Oncogene C (Gain): {oncogene_C_gain} copies * {oncogene_C_weight}/copy = {score_C}")
print(f"Tumor Suppressor D (Repressor-induced loss): {tsg_D_functional_loss} copies * {tsg_D_weight}/copy = {score_D}")
print("-" * 20)
print(f"Final Equation: {score_A} + {score_B} + {score_C} + {score_D} = {total_score}")
print(f"Total Clonal Expansion Score: {total_score}")