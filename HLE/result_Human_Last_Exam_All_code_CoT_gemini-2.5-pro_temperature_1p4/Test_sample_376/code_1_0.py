# Define the copy number variations (CNVs)
# Positive for gain, negative for loss
cnv_chr1 = 3  # Gain of 3 copies
cnv_chr2 = -2 # Loss of 2 copies
cnv_chr3 = 2  # Gain of 2 copies

# Define the weights for the impact on clonal expansion
weight_oncogene_A = 0.5   # per additional copy
weight_ts_B = -0.7        # per lost copy
weight_oncogene_C = 0.4   # per additional copy
# The following weights are for events that do not occur, but are defined for completeness.
# weight_ts_D = -0.6      # per lost copy
# weight_oncogene_E = 0.3 # per additional copy
# weight_ts_F = -0.5      # per lost copy

# --- Calculation ---

# Chromosome 1: Gain of 3 copies
# Oncogene A is gained. Its copy number increases by 3.
score_A = cnv_chr1 * weight_oncogene_A
# Tumor suppressor D is gained. The weight is for loss, so its contribution is 0.
score_D = 0

# Chromosome 2: Loss of 2 copies
# Tumor suppressor B is lost. Its copy number decreases by 2.
# The loss of a tumor suppressor promotes expansion, so we use the absolute value of the change.
score_B = abs(cnv_chr2) * abs(weight_ts_B)
# Oncogene E is lost. The weight is for gain, so its contribution is 0.
score_E = 0

# Chromosome 3: Gain of 2 copies
# Oncogene C is gained. Its copy number increases by 2.
score_C = cnv_chr3 * weight_oncogene_C
# Tumor suppressor F is gained. The weight is for loss, so its contribution is 0.
score_F = 0

# The repressor information is not used as no numerical weight is provided.

# Calculate the total clonal expansion score
total_score = score_A + score_B + score_C + score_D + score_E + score_F

# Print the final equation and the result
print("Clonal Expansion Score Calculation:")
print(f"{score_A} (from Oncogene A) + {score_B:.1f} (from Tumor Suppressor B) + {score_C:.1f} (from Oncogene C) = {total_score:.1f}")
