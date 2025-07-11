# Define data for each user connecting to B
# P1: 7 total, 6 trust, 1 distrust -> gives Trust edge to B
p1_total = 7
# P2: 6 total, 4 trust, 2 distrust -> gives Trust edge to B
p2_total = 6
# P3: 4 total, 2 trust, 2 distrust -> gives Trust edge to B
p3_total = 4
# N1: 6 total, 3 trust, 3 distrust -> gives Distrust edge to B
n1_total = 6
n1_trust = 3
n1_distrust = 3
# N2: 4 total, 1 trust, 3 distrust -> gives Distrust edge to B
n2_total = 4
n2_trust = 1
n2_distrust = 3

# --- Calculations ---

# Step 1: Calculate contribution from positive edges (Rule 1)
p1_score = 1 / (p1_total + 1)
p2_score = 1 / (p2_total + 1)
p3_score = 1 / (p3_total + 1)

# Step 2: Calculate contribution from negative edges (Rules 2 & 3)

# For N1
# Rule 2: Negative edge with mixed ratings
# Rule 3 does not apply as n1_distrust is not > n1_trust
n1_score = -1 / (n1_total + 1) * (n1_trust / (n1_total + 1))

# For N2
# Rule 2: Negative edge with mixed ratings
# Rule 3 applies as n2_distrust > n2_trust, so multiply by 1.5
n2_score = 1.5 * (-1 / (n2_total + 1) * (n2_trust / (n2_total + 1)))

# Step 3: Sum all contributions
total_score = p1_score + p2_score + p3_score + n1_score + n2_score

# --- Output the results ---

print("Calculating B's importance score (FAT):")
print("\nContribution from Trust edges:")
print(f"P1: 1 / ({p1_total} + 1) = {p1_score:.4f}")
print(f"P2: 1 / ({p2_total} + 1) = {p2_score:.4f}")
print(f"P3: 1 / ({p3_total} + 1) = {p3_score:.4f}")

print("\nContribution from Distrust edges:")
print(f"N1: -1 / ({n1_total} + 1) * ({n1_trust} / ({n1_total} + 1)) = {n1_score:.4f}")
print(f"N2: 1.5 * [-1 / ({n2_total} + 1) * ({n2_trust} / ({n2_total} + 1))] = {n2_score:.4f}")

print("\nFinal Equation:")
print(f"Score = {p1_score:.4f} + {p2_score:.4f} + {p3_score:.4f} + ({n1_score:.4f}) + ({n2_score:.4f})")

print(f"\nTotal Score = {total_score:.4f}")
print(f"The closest answer is 0.35.")
