import math

# Define the data for each user connected to B.
# Positive connections (Trust Incoming)
p1_data = {'total': 7}
p2_data = {'total': 6}
p3_data = {'total': 4}

# Negative connections (Distrust Incoming)
n1_data = {'total': 6, 'trust': 3, 'distrust': 3}
n2_data = {'total': 4, 'trust': 1, 'distrust': 3}

# --- Calculation based on FAT rules ---

# Rule 1: Positive edge contributes 1/(total_relationships + 1)
p1_contrib = 1 / (p1_data['total'] + 1)
p2_contrib = 1 / (p2_data['total'] + 1)
p3_contrib = 1 / (p3_data['total'] + 1)

# Rule 2 & 3 for negative edges
# N1's contribution
n1_contrib = -1 / (n1_data['total'] + 1) * (n1_data['trust'] / n1_data['total'])
# For N1, distrust (3) is not > trust (3), so no 1.5x multiplier.

# N2's contribution
n2_contrib = -1 / (n2_data['total'] + 1) * (n2_data['trust'] / n2_data['total'])
# For N2, distrust (3) > trust (1), so apply 1.5x multiplier from Rule 3.
n2_contrib *= 1.5

# Sum all contributions for the final score
total_score = p1_contrib + p2_contrib + p3_contrib + n1_contrib + n2_contrib

# --- Output the results ---
print("Calculating B's Fringe of Absolute Trust (FAT) score:")
print("-" * 50)
print(f"P1 (Positive): 1 / ({p1_data['total']} + 1) = {p1_contrib:.4f}")
print(f"P2 (Positive): 1 / ({p2_data['total']} + 1) = {p2_contrib:.4f}")
print(f"P3 (Positive): 1 / ({p3_data['total']} + 1) = {p3_contrib:.4f}")
print(f"N1 (Negative): -1 / ({n1_data['total']} + 1) * ({n1_data['trust']}/{n1_data['total']}) = {n1_contrib:.4f}")
print(f"N2 (Negative): [-1 / ({n2_data['total']} + 1) * ({n2_data['trust']}/{n2_data['total']})] * 1.5 = {n2_contrib:.4f}")
print("-" * 50)
print("Final Equation:")
print(f"{p1_contrib:.4f} + {p2_contrib:.4f} + {p3_contrib:.4f} + ({n1_contrib:.4f}) + ({n2_contrib:.4f}) = {total_score:.4f}")
print(f"\nThe calculated score is {total_score:.4f}, which is closest to 0.35.")