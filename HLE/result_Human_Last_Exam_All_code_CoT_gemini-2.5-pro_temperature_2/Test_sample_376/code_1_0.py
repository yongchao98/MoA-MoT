# Data from the problem description
# CNV data
gain_chr1 = 3
loss_chr2 = 2
gain_chr3 = 2

# Gene impact weights
w_oncogene_a = 0.5 # per additional copy
w_tsg_b = -0.7      # per lost copy
w_oncogene_c = 0.4 # per additional copy
# Note: Other genes (TSG D, Oncogene E, TSG F) have weights for changes
# (loss, gain, loss respectively) that did not occur, so their impact is 0.

# --- Calculation ---
# Chromosome 1: Gain of 3 affects Oncogene A (gain) and TSG D (gain)
# Oncogene A gain is pro-expansion. TSG D gain is not a loss, so its impact is 0.
score_oncogene_a = gain_chr1 * w_oncogene_a

# Chromosome 2: Loss of 2 affects TSG B (loss) and Oncogene E (loss)
# TSG B loss is pro-expansion. Oncogene E loss is not a gain, so its impact is 0.
# The impact score for a loss of a tumor suppressor should be positive.
score_tsg_b = loss_chr2 * abs(w_tsg_b)

# Chromosome 3: Gain of 2 affects Oncogene C (gain) and TSG F (gain)
# Oncogene C gain is pro-expansion. TSG F gain is not a loss, so its impact is 0.
score_oncogene_c = gain_chr3 * w_oncogene_c

# Total clonal expansion score
total_score = score_oncogene_a + score_tsg_b + score_oncogene_c

# --- Output ---
print("The clonal expansion score is calculated by summing the impacts of each relevant genetic event.")
print("Impacts are calculated as (number of copies changed) * (weight per copy).")
print("Only changes that have a defined weight contribute to the score.")
print("\nContributing factors:")
print(f"- Oncogene A gain: {gain_chr1} copies * {w_oncogene_a} weight = {score_oncogene_a}")
print(f"- Tumor Suppressor B loss: {loss_chr2} copies * {abs(w_tsg_b)} weight = {score_tsg_b:.1f}")
print(f"- Oncogene C gain: {gain_chr3} copies * {w_oncogene_c} weight = {score_oncogene_c:.1f}")
print("\nFinal Equation:")
print(f"({gain_chr1} * {w_oncogene_a}) + ({loss_chr2} * {abs(w_tsg_b)}) + ({gain_chr3} * {w_oncogene_c}) = {total_score}")

print(f"\nThe total clonal expansion score is: {total_score}")
<<<3.7>>>