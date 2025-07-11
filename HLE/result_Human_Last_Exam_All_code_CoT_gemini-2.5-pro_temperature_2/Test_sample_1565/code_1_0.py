# User data based on the problem description
# P_users have positive (trust) edges to B
# N_users have negative (distrust) edges to B

# P1: 7 total ratings, 6 trust, 1 distrust
p1_total = 7

# P2: 6 total ratings, 4 trust, 2 distrust
p2_total = 6

# P3: 4 total ratings, 2 trust, 2 distrust
p3_total = 4

# N1: 6 total ratings, 3 trust, 3 distrust
n1_total = 6
n1_trust = 3
n1_distrust = 3

# N2: 4 total ratings, 1 trust, 3 distrust
n2_total = 4
n2_trust = 1
n2_distrust = 3

# --- Calculation Step-by-Step ---

# Rule 1: Positive edge contributes 1 / (total_relationships + 1)
p1_contrib = 1 / (p1_total + 1)
p2_contrib = 1 / (p2_total + 1)
p3_contrib = 1 / (p3_total + 1)

# Rule 2: Negative edge with mixed ratings: -1 / (total_relationships + 1) * (trust_ratings / total)
n1_contrib = -1 / (n1_total + 1) * (n1_trust / n1_total)

# Rule 3: Users with more distrust than trust get 1.5x negative weight
n2_contrib_base = -1 / (n2_total + 1) * (n2_trust / n2_total)
n2_contrib = n2_contrib_base
if n2_distrust > n2_trust:
    n2_contrib *= 1.5

# --- Summing the scores ---
total_score = p1_contrib + p2_contrib + p3_contrib + n1_contrib + n2_contrib

# --- Output the results ---
print("Calculating B's importance score (FAT):")
print("\nContributions from each user:")
print(f"P1 (positive): 1 / ({p1_total} + 1) = {p1_contrib:.4f}")
print(f"P2 (positive): 1 / ({p2_total} + 1) = {p2_contrib:.4f}")
print(f"P3 (positive): 1 / ({p3_total} + 1) = {p3_contrib:.4f}")
print(f"N1 (negative): -1 / ({n1_total} + 1) * ({n1_trust} / {n1_total}) = {n1_contrib:.4f}")
if n2_distrust > n2_trust:
    print(f"N2 (negative, with 1.5x weight): (-1 / ({n2_total} + 1) * ({n2_trust} / {n2_total})) * 1.5 = {n2_contrib:.4f}")
else:
    print(f"N2 (negative): -1 / ({n2_total} + 1) * ({n2_trust} / {n2_total}) = {n2_contrib:.4f}")

print("\nFinal Equation:")
# The formatting below shows each calculated number in the final equation as requested
final_equation = f"{p1_contrib:.4f} + {p2_contrib:.4f} + {p3_contrib:.4f} + ({n1_contrib:.4f}) + ({n2_contrib:.4f}) = {total_score:.4f}"
print(final_equation)

print(f"\nThe calculated score is {total_score:.4f}, which is closest to 0.35.")
print("A")