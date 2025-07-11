# Define the given copy number changes and weights
gain_chr1 = 3
loss_chr2 = -2
gain_chr3 = 2

weight_oncogene_A = 0.5
weight_ts_B = -0.7
weight_oncogene_C = 0.4
weight_ts_D = -0.6
# Oncogene E and Tumor Suppressor F weights are not needed as their CNV events
# (loss of an oncogene, gain of a tumor suppressor) do not have corresponding scoring rules.

# --- Calculate the score contribution from each relevant event ---

# 1. Oncogene A (Chr 1): Score from gain of 3 copies
score_oncogene_A = gain_chr1 * weight_oncogene_A

# 2. Tumor Suppressor B (Chr 2): Score from loss of 2 copies
# A loss of 2 copies (-2) multiplied by a negative weight (-0.7) results in a positive score,
# correctly reflecting that losing a tumor suppressor promotes expansion.
score_ts_B = loss_chr2 * weight_ts_B

# 3. Oncogene C (Chr 3): Score from gain of 2 copies
score_oncogene_C = gain_chr3 * weight_oncogene_C

# 4. Tumor Suppressor D (Chr 1): Score from repressor action
# The gene has a copy number gain, which does not contribute directly to the score.
# However, the repressor causes a functional loss of the normal 2 copies.
functional_loss_D = -2
score_ts_D = functional_loss_D * weight_ts_D

# --- Calculate the total score and print the final equation ---

total_score = score_oncogene_A + score_ts_B + score_oncogene_C + score_ts_D

print("Clonal Expansion Score Calculation:")
# The following line prints each individual score that contributes to the final sum.
print(f"{score_oncogene_A} (from Oncogene A) + {score_ts_B} (from Tumor Suppressor B) + {score_oncogene_C} (from Oncogene C) + {score_ts_D} (from Tumor Suppressor D) = {total_score}")