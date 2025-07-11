# Define the copy number changes
chr1_gain = 3
chr2_loss = 2
chr3_gain = 2

# Define the weights for the impact of each gene
weight_oncogene_A = 0.5  # per additional copy
weight_tsg_B = -0.7       # per lost copy
weight_oncogene_C = 0.4  # per additional copy
# Note: Other weights are not applicable to the given CNVs (e.g., weight for lost copy of TSG D when it's a gain)

# --- Calculate the contribution of each relevant event ---

# Chromosome 1: Gain of 3 copies affects Oncogene A
score_oncogene_A = chr1_gain * weight_oncogene_A

# Chromosome 2: Loss of 2 copies affects Tumor Suppressor B.
# A loss of a tumor suppressor increases the clonal expansion score, so we use the absolute value of the weight.
score_tsg_B = chr2_loss * abs(weight_tsg_B)

# Chromosome 3: Gain of 2 copies affects Oncogene C
score_oncogene_C = chr3_gain * weight_oncogene_C

# --- Calculate the total clonal expansion score ---
total_score = score_oncogene_A + score_tsg_B + score_oncogene_C

# --- Print the final equation and result ---
# The format requires printing each number in the final equation.
print("Clonal Expansion Score Calculation:")
print(f"{score_oncogene_A} (from Oncogene A) + {score_tsg_B:.1f} (from Tumor Suppressor B) + {score_oncogene_C} (from Oncogene C) = {total_score:.1f}")