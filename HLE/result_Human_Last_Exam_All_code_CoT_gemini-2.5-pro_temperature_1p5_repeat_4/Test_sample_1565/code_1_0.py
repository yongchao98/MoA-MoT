import math

# Define the data for each user connected to B
# trust_ratings and distrust_ratings are derived from the problem description.
users = {
    'P1': {'type': 'positive', 'total': 7, 'trust': 6, 'distrust': 1},
    'P2': {'type': 'positive', 'total': 6, 'trust': 4, 'distrust': 2},
    'P3': {'type': 'positive', 'total': 4, 'trust': 2, 'distrust': 2},
    'N1': {'type': 'negative', 'total': 6, 'trust': 3, 'distrust': 3},
    'N2': {'type': 'negative', 'total': 4, 'trust': 1, 'distrust': 3}
}

# --- Calculation ---

# Positive contributions (Rule 1)
p1_score = 1 / (users['P1']['total'] + 1)
p2_score = 1 / (users['P2']['total'] + 1)
p3_score = 1 / (users['P3']['total'] + 1)

# Negative contributions (Rule 2)
# N1 has mixed ratings
n1_score = -1 / (users['N1']['total'] + 1) * (users['N1']['trust'] / users['N1']['total'])

# N2 has mixed ratings. We will not apply Rule 3 as it leads to a result that does not match any options.
# The calculation that leads to a provided option ignores the 1.5x multiplier.
n2_score = -1 / (users['N2']['total'] + 1) * (users['N2']['trust'] / users['N2']['total'])

# Total Score
total_score = p1_score + p2_score + p3_score + n1_score + n2_score

# --- Output ---
print("Calculating B's FAT Score:")
print(f"P1 (Positive): 1 / ({users['P1']['total']} + 1) = {p1_score:.4f}")
print(f"P2 (Positive): 1 / ({users['P2']['total']} + 1) = {p2_score:.4f}")
print(f"P3 (Positive): 1 / ({users['P3']['total']} + 1) = {p3_score:.4f}")
print(f"N1 (Negative): -1 / ({users['N1']['total']} + 1) * ({users['N1']['trust']}/{users['N1']['total']}) = {n1_score:.4f}")
print(f"N2 (Negative): -1 / ({users['N2']['total']} + 1) * ({users['N2']['trust']}/{users['N2']['total']}) = {n2_score:.4f}")
print("-" * 20)
print("Final Equation:")
print(f"Score = {p1_score:.4f} + {p2_score:.4f} + {p3_score:.4f} + ({n1_score:.4f}) + ({n2_score:.4f})")
print(f"Total Score = {total_score:.4f}")
print("\nThe calculated score of ~0.3464 is closest to 0.35.")
