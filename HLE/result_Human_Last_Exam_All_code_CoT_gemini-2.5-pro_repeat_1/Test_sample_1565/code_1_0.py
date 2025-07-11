import math

# Step 1: Define the data for each user connecting to B.
# P1, P2, P3 have positive (trust) edges to B.
# N1, N2 have negative (distrust) edges to B.
users = {
    'P1': {'total_ratings': 7},
    'P2': {'total_ratings': 6},
    'P3': {'total_ratings': 4},
    'N1': {'total_ratings': 6, 'trust_ratings': 3, 'distrust_ratings': 3},
    'N2': {'total_ratings': 4, 'trust_ratings': 1, 'distrust_ratings': 3}
}

# Step 2: Calculate the contribution of each user.

# For positive edges (P1, P2, P3), apply Rule 1.
p1_score = 1 / (users['P1']['total_ratings'] + 1)
p2_score = 1 / (users['P2']['total_ratings'] + 1)
p3_score = 1 / (users['P3']['total_ratings'] + 1)

# For negative edge from N1, apply Rule 2.
# N1 has mixed ratings (3 trust, 3 distrust).
# Distrust is not greater than trust, so Rule 3 does not apply.
n1_data = users['N1']
n1_score = -1 / (n1_data['total_ratings'] + 1) * (n1_data['trust_ratings'] / n1_data['total_ratings'])

# For negative edge from N2, apply Rule 2 and Rule 3.
# N2 has mixed ratings (1 trust, 3 distrust).
# Distrust (3) is greater than trust (1), so Rule 3 (1.5x weight) applies.
n2_data = users['N2']
n2_score = 1.5 * (-1 / (n2_data['total_ratings'] + 1) * (n2_data['trust_ratings'] / n2_data['total_ratings']))

# Step 3: Sum the contributions for the final score.
total_score = p1_score + p2_score + p3_score + n1_score + n2_score

# Output the final equation with each number and the result.
# The problem asks to show each number in the final equation.
print(f"{p1_score:.4f} + {p2_score:.4f} + {p3_score:.4f} + ({n1_score:.4f}) + ({n2_score:.4f}) = {total_score:.4f}")

# The calculated score is ~0.3214, which is closest to option A (0.35).
print("A")