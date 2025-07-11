# Define the copy number changes and weights for relevant genes.
# A "gain" of an oncogene contributes positively to the score.
# A "loss" of a tumor suppressor also contributes positively, as it aids cancer progression.

# Chromosome 1: Gain of 3 copies affecting Oncogene A
gain_oncogene_A = 3
weight_oncogene_A = 0.5

# Chromosome 2: Loss of 2 copies affecting Tumor Suppressor B
loss_ts_B = 2
# The weight is given as -0.7 per lost copy. We interpret this as a positive
# contribution to clonal expansion, so we use 0.7.
weight_ts_B = 0.7

# Chromosome 3: Gain of 2 copies affecting Oncogene C
gain_oncogene_C = 2
weight_oncogene_C = 0.4

# Calculate the score contribution from each event
score_A = gain_oncogene_A * weight_oncogene_A
score_B = loss_ts_B * weight_ts_B
score_C = gain_oncogene_C * weight_oncogene_C

# The other genes (D, E, F) do not contribute to the score because the type of
# change they underwent (gain of TS, loss of oncogene) has no associated weight.

# Calculate the final clonal expansion score
total_score = score_A + score_B + score_C

# Print the final equation, showing each number involved in the calculation
print(f"Clonal Expansion Score = (Contribution from Oncogene A) + (Contribution from Tumor Suppressor B) + (Contribution from Oncogene C)")
print(f"Clonal Expansion Score = ({gain_oncogene_A} * {weight_oncogene_A}) + ({loss_ts_B} * {weight_ts_B}) + ({gain_oncogene_C} * {weight_oncogene_C})")
print(f"Clonal Expansion Score = {score_A} + {score_B} + {score_C} = {total_score}")
<<<3.7>>>