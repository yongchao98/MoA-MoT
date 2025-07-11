# Define the copy number changes
chr1_gain = 3
chr2_loss = 2
chr3_gain = 2

# Define the weights for the impact of each gene
weight_onc_A = 0.5  # per additional copy
weight_ts_B = -0.7  # per lost copy
weight_onc_C = 0.4  # per additional copy

# Calculate the score from each relevant event
# Oncogene A gain on Chromosome 1
score_onc_A = chr1_gain * weight_onc_A

# Tumor suppressor B loss on Chromosome 2
# The loss of a suppressor promotes expansion, so we use the absolute value of the weight
score_ts_B = chr2_loss * abs(weight_ts_B)

# Oncogene C gain on Chromosome 3
score_onc_C = chr3_gain * weight_onc_C

# Other genes (TS D, Onc E, TS F) do not contribute to the score based on the provided weights and events
# (e.g., TS D is gained but weight is for loss, Onc E is lost but weight is for gain).

# Calculate the total clonal expansion score
total_score = score_onc_A + score_ts_B + score_onc_C

# Print the final equation and the result
print("The clonal expansion score is calculated by summing the impact of each relevant genetic event.")
print(f"Final Score = (Score from Oncogene A) + (Score from Tumor Suppressor B) + (Score from Oncogene C)")
print(f"Final Equation: {score_onc_A} + {score_ts_B} + {score_onc_C} = {total_score}")
