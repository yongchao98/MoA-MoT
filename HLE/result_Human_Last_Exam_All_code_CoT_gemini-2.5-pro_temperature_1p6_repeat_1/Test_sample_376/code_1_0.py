import math

# Plan:
# 1. Define variables for copy number changes and weights.
# 2. Interpret weights: A positive contribution to clonal expansion score occurs
#    from oncogene gains and tumor suppressor losses. The magnitude of the weight
#    will be used for the calculation.
# 3. Calculate the contribution for each gene only if the CNV event (gain/loss)
#    matches the condition for its weight.
# 4. Sum the contributions and print the detailed calculation.

# Chromosome Copy Number Variations (CNVs)
chr1_cnv = 3  # gain of 3
chr2_cnv = -2 # loss of 2
chr3_cnv = 2  # gain of 2

# Weights for impact on clonal expansion
weight_oncogene_A = 0.5  # per additional copy
weight_tsg_B = -0.7       # per lost copy
weight_oncogene_C = 0.4  # per additional copy
weight_tsg_D = -0.6       # per lost copy
weight_oncogene_E = 0.3  # per additional copy
weight_tsg_F = -0.5       # per lost copy

# Initialize scores for each gene
score_A = 0
score_B = 0
score_C = 0
score_D = 0
score_E = 0
score_F = 0

print("Calculating the clonal expansion score based on genetic variations:")

# --- Gene-by-Gene Calculation ---

# Oncogene A on Chromosome 1 (gain of 3 copies)
# The weight applies to additional copies, and a gain occurred.
if chr1_cnv > 0:
    score_A = chr1_cnv * weight_oncogene_A
    print(f"- Contribution from Oncogene A (Gain of {chr1_cnv} copies * {weight_oncogene_A}/copy) = {score_A}")

# Tumor Suppressor D on Chromosome 1 (gain of 3 copies)
# The weight is for lost copies, but a gain occurred, so no contribution.
# No print output needed as contribution is 0.

# Tumor Suppressor B on Chromosome 2 (loss of 2 copies)
# The weight applies to lost copies, and a loss occurred. This promotes expansion.
if chr2_cnv < 0:
    num_lost = abs(chr2_cnv)
    # Use the absolute value of the weight for a positive score contribution
    score_B = num_lost * abs(weight_tsg_B)
    print(f"- Contribution from Tumor Suppressor B (Loss of {num_lost} copies * {abs(weight_tsg_B)}/copy) = {score_B:.1f}")

# Oncogene E on Chromosome 2 (loss of 2 copies)
# The weight is for additional copies, but a loss occurred, so no contribution.
# No print output needed as contribution is 0.

# Oncogene C on Chromosome 3 (gain of 2 copies)
# The weight applies to additional copies, and a gain occurred.
if chr3_cnv > 0:
    score_C = chr3_cnv * weight_oncogene_C
    print(f"- Contribution from Oncogene C (Gain of {chr3_cnv} copies * {weight_oncogene_C}/copy) = {score_C:.1f}")

# Tumor Suppressor F on Chromosome 3 (gain of 2 copies)
# The weight is for lost copies, but a gain occurred, so no contribution.
# No print output needed as contribution is 0.

# Sum the scores from all relevant events
total_score = score_A + score_B + score_C + score_D + score_E + score_F

print("- Other genes do not contribute as the CNV event did not match their scoring condition.")
print("\nFinal Equation:")

# Get the list of non-zero scores for the final equation string
active_scores = [s for s in [score_A, score_B, score_C] if s != 0]

# Format the numbers to avoid floating point inaccuracies in the string
formatted_scores = [f"{s:.1f}" if not s.is_integer() else str(int(s)) for s in active_scores]

equation_str = " + ".join(formatted_scores)

# Print the final equation with each contributing number
print(f"{equation_str} = {total_score:.1f}")

print(f"\nThe total clonal expansion score is: {total_score:.1f}")