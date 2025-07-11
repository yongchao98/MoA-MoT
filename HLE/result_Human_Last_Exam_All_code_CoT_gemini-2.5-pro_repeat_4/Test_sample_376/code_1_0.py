# Plan:
# 1. Define variables for the copy number changes and weights based on the problem description.
#    - Gains are positive numbers, losses are negative numbers.
# 2. Calculate the score contribution for each gene. A contribution is only calculated if the type of
#    copy number change (gain/loss) matches the condition provided for the weight.
# 3. Sum the individual contributions to get the final clonal expansion score.
# 4. Print the final calculation, showing each term in the equation.

# Chromosome 1: Gain of 3 copies
c1_change = 3
# Chromosome 2: Loss of 2 copies
c2_change = -2
# Chromosome 3: Gain of 2 copies
c3_change = 2

# Weights for impact on clonal expansion
w_onc_A = 0.5   # per additional copy
w_tsg_B = -0.7  # per lost copy
w_onc_C = 0.4   # per additional copy
# Note: Other weights are not used as the conditions are not met.

# Calculate contribution from each relevant gene
# Oncogene A: Gain of 3, weight is for gain. Contribution is calculated.
score_A = c1_change * w_onc_A

# Tumor suppressor B: Loss of 2, weight is for loss. Contribution is calculated.
score_B = c2_change * w_tsg_B

# Oncogene C: Gain of 2, weight is for gain. Contribution is calculated.
score_C = c3_change * w_onc_C

# The other genes do not contribute to the score:
# - Tumor suppressor D: Gained, but weight is for loss.
# - Oncogene E: Lost, but weight is for gain.
# - Tumor suppressor F: Gained, but weight is for loss.

# Calculate the total score
total_score = score_A + score_B + score_C

# Print the final equation with all numbers
print("The clonal expansion score is calculated by summing the impacts of each relevant genetic change.")
print("The calculation is as follows:")
print(f"Score = (Oncogene A contribution) + (Tumor Suppressor B contribution) + (Oncogene C contribution)")
print(f"Score = ({c1_change} * {w_onc_A}) + ({c2_change} * {w_tsg_B}) + ({c3_change} * {w_onc_C})")
print(f"Score = {score_A} + {score_B} + {score_C}")
print(f"Total Clonal Expansion Score = {total_score}")
