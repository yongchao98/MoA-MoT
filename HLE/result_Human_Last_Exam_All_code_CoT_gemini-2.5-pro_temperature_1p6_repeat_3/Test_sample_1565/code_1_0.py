# User data based on the problem description
# P for positive (trust) connections, N for negative (distrust)
# Data: (total_ratings, trust_ratings, distrust_ratings)
p1_data = {'total': 7}
p2_data = {'total': 6}
p3_data = {'total': 4}
n1_data = {'total': 6, 'trust': 3, 'distrust': 3}
n2_data = {'total': 4, 'trust': 1, 'distrust': 3}

# --- Calculation ---

# Rule 1: Positive edge contributes 1 / (total_relationships + 1)
p1_contrib = 1 / (p1_data['total'] + 1)
p2_contrib = 1 / (p2_data['total'] + 1)
p3_contrib = 1 / (p3_data['total'] + 1)

# Rule 2: Negative edge with mixed ratings: -1 / (total + 1) * (trust / total)
n1_contrib = -1 / (n1_data['total'] + 1) * (n1_data['trust'] / n1_data['total'])

# For N2, check Rule 3: Users with more distrust than trust get 1.5x negative weight
n2_contrib_base = -1 / (n2_data['total'] + 1) * (n2_data['trust'] / n2_data['total'])
n2_multiplier = 1.0
if n2_data['distrust'] > n2_data['trust']:
    n2_multiplier = 1.5
n2_contrib = n2_contrib_base * n2_multiplier

# --- Summing the scores ---
total_score = p1_contrib + p2_contrib + p3_contrib + n1_contrib + n2_contrib

# --- Output the results ---
# Displaying each part of the equation as requested.
p1_str = f"(1 / ({p1_data['total']} + 1))"
p2_str = f"(1 / ({p2_data['total']} + 1))"
p3_str = f"(1 / ({p3_data['total']} + 1))"
n1_str = f"(-1 / ({n1_data['total']} + 1) * ({n1_data['trust']} / {n1_data['total']}))"
n2_str = f"((-1 / ({n2_data['total']} + 1) * ({n2_data['trust']} / {n2_data['total']})) * {n2_multiplier})"

print("B's FAT Score Calculation:")
print(f"Contribution from P1: {p1_str} = {p1_contrib:.4f}")
print(f"Contribution from P2: {p2_str} = {p2_contrib:.4f}")
print(f"Contribution from P3: {p3_str} = {p3_contrib:.4f}")
print(f"Contribution from N1: {n1_str} = {n1_contrib:.4f}")
print(f"Contribution from N2: {n2_str} = {n2_contrib:.4f}")
print("\nFinal Equation:")
# We re-format the numbers to match the step-by-step part
# for clarity in the final equation printout
p1_eq = f"(1 / 8)"
p2_eq = f"(1 / 7)"
p3_eq = f"(1 / 5)"
n1_eq = f"(-1 / 7 * 3 / 6)"
n2_eq = f"((-1 / 5 * 1 / 4) * 1.5)"

print(f"Score = {p1_eq} + {p2_eq} + {p3_eq} + {n1_eq} + {n2_eq}")

print(f"\nTotal Score = {p1_contrib:.4f} + {p2_contrib:.4f} + {p3_contrib:.4f} + {n1_contrib:.4f} + {n2_contrib:.4f} = {total_score:.4f}")

# The calculated score is ~0.3214. The closest answer choice is 0.35.
print("\nThe calculated score is approximately 0.3214.")
print("The closest option is A) 0.35.")
