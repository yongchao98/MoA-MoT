# Step 1: Define user data
# (total_ratings, trust_ratings)
p1_data = (7, 6)
p2_data = (6, 4)
p3_data = (4, 2)
n1_data = (6, 3)
n2_data = (4, 1)

# Step 2: Calculate contribution for each user
# Positive edges (Rule 1)
p1_contrib = 1 / (p1_data[0] + 1)
p2_contrib = 1 / (p2_data[0] + 1)
p3_contrib = 1 / (p3_data[0] + 1)

# Negative edge N1 (Rule 2)
n1_contrib = -1 / (n1_data[0] + 1) * (n1_data[1] / n1_data[0])

# Negative edge N2 (Rules 2 and 3)
# N2 has 1 trust and 3 distrust, so distrust > trust. Apply 1.5x multiplier.
n2_contrib = (-1 / (n2_data[0] + 1) * (n2_data[1] / n2_data[0])) * 1.5

# Step 3: Sum the contributions
total_score = p1_contrib + p2_contrib + p3_contrib + n1_contrib + n2_contrib

# Output the final equation and result
print("B's importance score calculation:")
print(f"({p1_contrib:.3f} from P1) + ({p2_contrib:.3f} from P2) + ({p3_contrib:.3f} from P3) + ({n1_contrib:.3f} from N1) + ({n2_contrib:.3f} from N2)")
print(f"= (1/8) + (1/7) + (1/5) - (1/14) - (3/40)")
print(f"= {p1_contrib} + {p2_contrib} + {p3_contrib} + {n1_contrib} + {n2_contrib}")
print(f"Total Score = {total_score}")
print(f"The calculated score is approximately {total_score:.2f}. The closest option is 0.35.")
print("\nA")