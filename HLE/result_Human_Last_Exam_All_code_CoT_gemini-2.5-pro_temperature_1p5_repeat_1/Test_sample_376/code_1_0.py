# --- Define the copy number changes and weights ---

# Chromosome 1: Gain of 3 copies
chr1_gain = 3
oncogene_A_weight = 0.5

# Chromosome 2: Loss of 2 copies
chr2_loss = 2
# The loss of a tumor suppressor promotes expansion, so we use the magnitude of the weight.
tsg_B_weight = 0.7

# Chromosome 3: Gain of 2 copies
chr3_gain = 2
oncogene_C_weight = 0.4

# --- Calculate the contribution of each relevant event to the score ---

# Contribution from Oncogene A gain on Chromosome 1
score_oncogene_A = chr1_gain * oncogene_A_weight

# Contribution from Tumor Suppressor B loss on Chromosome 2
score_tsg_B = chr2_loss * tsg_B_weight

# Contribution from Oncogene C gain on Chromosome 3
score_oncogene_C = chr3_gain * oncogene_C_weight

# Other genes (Tumor suppressor D, Oncogene E, Tumor suppressor F) do not contribute
# as the specific events affecting them (gain of TSG, loss of Oncogene) have no assigned weights.

# --- Calculate the total clonal expansion score ---
total_score = score_oncogene_A + score_tsg_B + score_oncogene_C

# --- Print the detailed calculation and the final result ---
print("The total clonal expansion score is calculated by summing the impact of each relevant genetic event.")
print("The final equation is derived as follows:")
print(f"({chr1_gain} * {oncogene_A_weight}) + ({chr2_loss} * {tsg_B_weight}) + ({chr3_gain} * {oncogene_C_weight}) = {total_score}")
