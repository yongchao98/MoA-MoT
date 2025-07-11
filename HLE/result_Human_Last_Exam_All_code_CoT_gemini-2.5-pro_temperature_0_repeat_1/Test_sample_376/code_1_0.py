# Define the copy number changes and weights based on the problem description.
# For tumor suppressors, a loss promotes clonal expansion, so we use a positive weight.
oncogene_A_gain = 3
weight_A = 0.5

tsg_B_loss = 2
weight_B = 0.7

oncogene_C_gain = 2
weight_C = 0.4

# Calculate the score contribution from each relevant gene
score_A = oncogene_A_gain * weight_A
score_B = tsg_B_loss * weight_B
score_C = oncogene_C_gain * weight_C

# Calculate the total clonal expansion score
total_score = score_A + score_B + score_C

# Print the final equation with each number and the total score
print("The clonal expansion score is calculated as follows:")
print(f"({oncogene_A_gain} * {weight_A}) + ({tsg_B_loss} * {weight_B}) + ({oncogene_C_gain} * {weight_C}) = {total_score}")
print("\nBreaking it down by gene:")
print(f"Oncogene A contribution: {score_A}")
print(f"Tumor Suppressor B contribution: {score_B}")
print(f"Oncogene C contribution: {score_C}")
print("\nFinal Equation:")
print(f"{score_A} + {score_B} + {score_C} = {total_score}")