# Step 1: Define user data based on the problem description.
# The format is: User: (total_relationships, trust_ratings, distrust_ratings)
users = {
    "P1": (7, 4 + 2, 1),      # Trust Incoming
    "P2": (6, 3 + 1, 2),      # Trust Incoming
    "P3": (4, 2, 1 + 1),      # Trust Incoming
    "N1": (6, 3, 2 + 1),      # Distrust Incoming
    "N2": (4, 1, 3)           # Distrust Incoming
}

# Initialize B's score
b_fat_score = 0

# Step 2: Calculate contribution for positive edges
print("--- Calculating Contributions ---")

# Contribution from P1
p1_total = users["P1"][0]
p1_contrib = 1 / (p1_total + 1)
b_fat_score += p1_contrib
print(f"P1 (Trust) contribution: 1 / ({p1_total} + 1) = {p1_contrib:.4f}")

# Contribution from P2
p2_total = users["P2"][0]
p2_contrib = 1 / (p2_total + 1)
b_fat_score += p2_contrib
print(f"P2 (Trust) contribution: 1 / ({p2_total} + 1) = {p2_contrib:.4f}")

# Contribution from P3
p3_total = users["P3"][0]
p3_contrib = 1 / (p3_total + 1)
b_fat_score += p3_contrib
print(f"P3 (Trust) contribution: 1 / ({p3_total} + 1) = {p3_contrib:.4f}")

# Step 3: Calculate contribution for negative edges

# Contribution from N1
n1_total = users["N1"][0]
n1_trust = users["N1"][1]
n1_distrust = users["N1"][2]
# Rule 3 check: N1's distrust (3) is not greater than trust (3), so no multiplier.
n1_contrib = -1 / (n1_total + 1) * (n1_trust / n1_total)
b_fat_score += n1_contrib
print(f"N1 (Distrust) contribution: -1 / ({n1_total} + 1) * ({n1_trust} / {n1_total}) = {n1_contrib:.4f}")

# Contribution from N2
n2_total = users["N2"][0]
n2_trust = users["N2"][1]
n2_distrust = users["N2"][2]
# Rule 3 check: N2's distrust (3) is greater than trust (1), so apply 1.5x multiplier.
n2_multiplier = 1.5
n2_contrib = -1 / (n2_total + 1) * (n2_trust / n2_total) * n2_multiplier
b_fat_score += n2_contrib
print(f"N2 (Distrust) contribution: -1 / ({n2_total} + 1) * ({n2_trust} / {n2_total}) * {n2_multiplier} = {n2_contrib:.4f}")

# Step 4 & 5: Sum contributions and present the final equation and result.
print("\n--- Final Score Calculation ---")
print("B's Importance Score = (P1 score) + (P2 score) + (P3 score) + (N1 score) + (N2 score)")
print(f"Final Equation: {p1_contrib:.4f} + {p2_contrib:.4f} + {p3_contrib:.4f} + ({n1_contrib:.4f}) + ({n2_contrib:.4f})")
print(f"Final Score = {b_fat_score:.4f}")