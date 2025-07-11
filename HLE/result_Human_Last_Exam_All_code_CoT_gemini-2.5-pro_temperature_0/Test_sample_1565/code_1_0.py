# 1. Define the data for each user connecting to B based on the problem description.
# P1, P2, P3 have positive edges (Trust Incoming).
# N1, N2 have negative edges (Distrust Incoming).
# Data includes total relationships, and counts of trust/distrust ratings.
p1_total_rels = 7
p2_total_rels = 6
p3_total_rels = 4

n1_total_rels = 6
n1_trust_ratings = 3
n1_distrust_ratings = 3 # 2 distrust-peer + 1 distrust-junior

n2_total_rels = 4
n2_trust_ratings = 1
n2_distrust_ratings = 3

# 2. Calculate the contribution from each user based on the FAT rules.

# Positive contributions (Rule 1: 1 / (total_relationships + 1))
p1_contrib = 1 / (p1_total_rels + 1)
p2_contrib = 1 / (p2_total_rels + 1)
p3_contrib = 1 / (p3_total_rels + 1)

# Negative contributions (Rule 2 & 3)
# N1: Negative edge with mixed ratings.
# Rule 3 check: distrust (3) is not > trust (3), so no 1.5x multiplier.
n1_contrib = -1 / (n1_total_rels + 1) * (n1_trust_ratings / n1_total_rels)

# N2: Negative edge with mixed ratings.
# Rule 3 check: distrust (3) > trust (1), so apply 1.5x multiplier.
n2_contrib_base = -1 / (n2_total_rels + 1) * (n2_trust_ratings / n2_total_rels)
n2_contrib = n2_contrib_base * 1.5

# 3. Sum all contributions to get the final FAT score for B.
total_score = p1_contrib + p2_contrib + p3_contrib + n1_contrib + n2_contrib

# 4. Print the detailed equation and the final result.
# The equation shows each number used in the calculation as requested.
equation = (
    f"B's score = (P1) + (P2) + (P3) + (N1) + (N2)\n"
    f"B's score = [1/({p1_total_rels}+1)] + [1/({p2_total_rels}+1)] + [1/({p3_total_rels}+1)] + "
    f"[-1/({n1_total_rels}+1) * ({n1_trust_ratings}/{n1_total_rels})] + "
    f"[-1/({n2_total_rels}+1) * ({n2_trust_ratings}/{n2_total_rels}) * 1.5]\n"
    f"B's score = {p1_contrib:.4f} + {p2_contrib:.4f} + {p3_contrib:.4f} + {n1_contrib:.4f} + {n2_contrib:.4f}"
)

print(equation)
print(f"Final FAT Score for B = {total_score:.4f}")

# The calculated score 0.3214 is closest to option A) 0.35.
print("\nThe calculated score is approximately 0.3214. Comparing this to the options:")
print("A) 0.35 (Difference: {:.4f})".format(abs(total_score - 0.35)))
print("B) 0.45 (Difference: {:.4f})".format(abs(total_score - 0.45)))
print("C) 0.55 (Difference: {:.4f})".format(abs(total_score - 0.55)))
print("D) 0.64 (Difference: {:.4f})".format(abs(total_score - 0.64)))
print("\nOption A is the closest answer.")