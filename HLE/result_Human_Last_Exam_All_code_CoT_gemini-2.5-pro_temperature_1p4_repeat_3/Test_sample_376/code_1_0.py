# Define the weights for each gene's impact
weights = {
    'Oncogene_A': 0.5,       # per additional copy
    'Tumor_suppressor_B': -0.7, # per lost copy
    'Oncogene_C': 0.4,       # per additional copy
    'Tumor_suppressor_D': -0.6, # per lost copy
    'Oncogene_E': 0.3,       # per additional copy
    'Tumor_suppressor_F': -0.5  # per lost copy
}

# --- Chromosome 1 Calculation ---
# Gain of 3 copies affects Oncogene A. This promotes clonal expansion.
gain_chr1 = 3
score_onc_A = gain_chr1 * weights['Oncogene_A']

# A repressor of Tumor suppressor D is overexpressed. This is functionally
# equivalent to losing the two normal copies of the gene, promoting expansion.
repressor_functional_loss = 2
score_ts_D_repressor = repressor_functional_loss * abs(weights['Tumor_suppressor_D'])

# --- Chromosome 2 Calculation ---
# Loss of 2 copies affects Tumor suppressor B. This promotes clonal expansion.
loss_chr2 = 2
score_ts_B = loss_chr2 * abs(weights['Tumor_suppressor_B'])

# Loss of 2 copies also affects Oncogene E. This inhibits clonal expansion.
score_onc_E = loss_chr2 * (-weights['Oncogene_E'])

# The repressor has no additional effect, as Tumor suppressor B is already lost.

# --- Chromosome 3 Calculation ---
# Gain of 2 copies affects Oncogene C. This promotes clonal expansion.
gain_chr3 = 2
score_onc_C = gain_chr3 * weights['Oncogene_C']

# --- Final Calculation ---
total_score = score_onc_A + score_ts_D_repressor + score_ts_B + score_onc_E + score_onc_C

print("The clonal expansion score is calculated by summing the impacts of all genetic events.")
print("\nContribution from each event:")
print(f"- Gain of Oncogene A: {score_onc_A}")
print(f"- Repression of Tumor Suppressor D: {score_ts_D_repressor}")
print(f"- Loss of Tumor Suppressor B: {score_ts_B}")
print(f"- Loss of Oncogene E: {score_onc_E}")
print(f"- Gain of Oncogene C: {score_onc_C}")

print("\nFinal Equation:")
# The format is used to represent the subtraction of a negative value as a simple subtraction
print(f"{score_onc_A} + {score_ts_D_repressor} + {score_ts_B} + ({score_onc_E}) + {score_onc_C} = {total_score}")

print(f"\nTotal Clonal Expansion Score: {total_score}")
