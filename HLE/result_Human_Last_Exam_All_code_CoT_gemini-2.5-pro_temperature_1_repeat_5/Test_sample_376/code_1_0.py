# Define the copy number variations for each chromosome
chr1_gain = 3
chr2_loss = 2
chr3_gain = 2

# Define the weights for the impact of each gene on clonal expansion
weight_oncogene_A = 0.5  # per additional copy
weight_suppressor_B = 0.7  # per lost copy (using positive magnitude for score)
weight_oncogene_C = 0.4  # per additional copy

# Calculate the impact score for each relevant genetic event
# Oncogene A on Chr 1 (gain of 3 copies)
score_A = chr1_gain * weight_oncogene_A

# Tumor suppressor B on Chr 2 (loss of 2 copies)
score_B = chr2_loss * weight_suppressor_B

# Oncogene C on Chr 3 (gain of 2 copies)
score_C = chr3_gain * weight_oncogene_C

# Tumor suppressor D, Oncogene E, and Tumor suppressor F do not contribute
# because the type of copy number change does not match their weight's condition.
# The repressor protein has no numerical weight, so it is not included.

# Calculate the total clonal expansion score
total_score = score_A + score_B + score_C

# Print the final equation with each number, as requested
print("Clonal Expansion Score Calculation:")
print(f"{score_A} (from Oncogene A) + {score_B:.1f} (from Tumor Suppressor B) + {score_C:.1f} (from Oncogene C) = {total_score:.1f}")
