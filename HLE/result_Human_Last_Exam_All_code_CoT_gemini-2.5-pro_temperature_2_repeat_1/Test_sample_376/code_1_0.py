# --- Define weights and copy number changes ---

# Chromosome 1: Gain of 3 copies
oncogene_A_additional_copies = 3
oncogene_A_weight = 0.5

# Chromosome 2: Loss of 2 copies
ts_B_lost_copies = 2
ts_B_weight = -0.7

# Chromosome 3: Gain of 2 copies
oncogene_C_additional_copies = 2
oncogene_C_weight = 0.4

# Repressor Effect on Tumor Suppressor D (on Chr 1)
# This is equivalent to losing the 2 normal copies of the gene.
ts_D_inactivated_copies = 2
ts_D_weight = -0.6

# --- Calculate the score for each component ---

# Gain of an oncogene is pro-expansion (positive score)
score_oncogene_A = oncogene_A_additional_copies * oncogene_A_weight

# Loss of a tumor suppressor is pro-expansion (positive score)
# The impact is positive, so we use the negative of the given weight.
score_ts_B = ts_B_lost_copies * (-ts_B_weight)

# Gain of an oncogene is pro-expansion (positive score)
score_oncogene_C = oncogene_C_additional_copies * oncogene_C_weight

# Repression of a tumor suppressor is equivalent to its loss, so it is pro-expansion (positive score)
score_repressor_ts_D = ts_D_inactivated_copies * (-ts_D_weight)

# --- Calculate the total clonal expansion score ---
total_score = score_oncogene_A + score_ts_B + score_oncogene_C + score_repressor_ts_D

# --- Print the final equation and result ---
# The equation shows each component of the calculation.
print("Clonal Expansion Score Calculation:")
print(
    f"({oncogene_A_additional_copies} * {oncogene_A_weight}) [Oncogene A gain] + "
    f"({ts_B_lost_copies} * {-ts_B_weight}) [Tumor Suppressor B loss] + "
    f"({oncogene_C_additional_copies} * {oncogene_C_weight}) [Oncogene C gain] + "
    f"({ts_D_inactivated_copies} * {-ts_D_weight}) [Repressor on TS D] = {total_score}"
)
<<<4.9>>>