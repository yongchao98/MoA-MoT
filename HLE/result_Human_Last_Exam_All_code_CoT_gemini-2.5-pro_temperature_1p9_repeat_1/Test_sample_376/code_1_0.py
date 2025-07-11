# --- Setup: Define variables from the problem statement ---

# Chromosomal copy number variations
chr1_gain = 3
chr2_loss = 2
chr3_gain = 2

# Weights for clonal expansion impact
# Note: A positive score contribution promotes expansion.
# Oncogene weights (per additional copy)
weight_oncogene_A = 0.5
weight_oncogene_C = 0.4

# Tumor Suppressor Gene (TSG) weights (per lost copy)
# The score *increases* for each lost copy.
# E.g., for TSG B, a loss adds 0.7 to the score.
score_per_loss_tsg_B = 0.7
score_per_loss_tsg_D = 0.6

# Repressor effect on chromosomes 1 and 2 is equivalent to a 1-copy loss
# of the respective TSG.
repressor_copy_loss_effect = 1

# --- Calculation: Step-by-step scoring ---

# Initialize a list to hold the individual score components for the equation
score_components = []
equation_parts = []

# Chromosome 1: Gain of 3 copies
# Impact of Oncogene A gain
score_A = chr1_gain * weight_oncogene_A
score_components.append(score_A)
equation_parts.append(f"({chr1_gain} * {weight_oncogene_A})")

# Impact of repressor on Tumor Suppressor D (equivalent to 1-copy loss)
score_D_repressor = repressor_copy_loss_effect * score_per_loss_tsg_D
score_components.append(score_D_repressor)
equation_parts.append(f"({repressor_copy_loss_effect} * {score_per_loss_tsg_D})")


# Chromosome 2: Loss of 2 copies
# Impact of Tumor Suppressor B copy loss
score_B_loss = chr2_loss * score_per_loss_tsg_B
score_components.append(score_B_loss)
equation_parts.append(f"({chr2_loss} * {score_per_loss_tsg_B})")

# Impact of repressor on Tumor Suppressor B (equivalent to an additional 1-copy loss)
score_B_repressor = repressor_copy_loss_effect * score_per_loss_tsg_B
score_components.append(score_B_repressor)
equation_parts.append(f"({repressor_copy_loss_effect} * {score_per_loss_tsg_B})")


# Chromosome 3: Gain of 2 copies
# Impact of Oncogene C gain
score_C = chr3_gain * weight_oncogene_C
score_components.append(score_C)
equation_parts.append(f"({chr3_gain} * {weight_oncogene_C})")

# --- Final Summation and Output ---
total_score = sum(score_components)

print("The clonal expansion score is the sum of the impacts from each genetic and regulatory event:")
print("\nScore = (Oncogene A Gain) + (TSG D Repressor) + (TSG B Loss) + (TSG B Repressor) + (Oncogene C Gain)")
print(f"Score = {' + '.join(equation_parts)}")
print(f"Score = {' + '.join(map(str, score_components))}")
print(f"\nTotal Clonal Expansion Score = {total_score}")
