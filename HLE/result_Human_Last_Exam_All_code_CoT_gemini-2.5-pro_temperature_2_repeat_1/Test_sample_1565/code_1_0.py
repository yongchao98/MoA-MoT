import math

# In a professional review network, user B has following relationships from users P1, P2, P3, N1, and N2.
# This section defines the rating profiles for each of these users.

# Users with "Trust Incoming" relationships to B (positive edges)
p1 = {'name': 'P1', 'total': 7, 'trust': 4 + 2, 'distrust': 1}
p2 = {'name': 'P2', 'total': 6, 'trust': 3 + 1, 'distrust': 2}
p3 = {'name': 'P3', 'total': 4, 'trust': 2, 'distrust': 1 + 1}

# Users with "Distrust Incoming" relationships to B (negative edges)
n1 = {'name': 'N1', 'total': 6, 'trust': 3, 'distrust': 2 + 1}
n2 = {'name': 'N2', 'total': 4, 'trust': 1, 'distrust': 3}

# --- Calculation based on FAT measure rules ---

# Rule 1: Positive edge contributes 1/(total_relationships + 1)
# Applied to P1, P2, P3
p1_contrib = 1 / (p1['total'] + 1)
p2_contrib = 1 / (p2['total'] + 1)
p3_contrib = 1 / (p3['total'] + 1)

# Rule 2: Negative edge with mixed ratings: -1 /(total relationships+1) x (trust_ratings/total)
# Applied to N1, as its distrust ratings are not greater than its trust ratings.
n1_contrib = -1 / (n1['total'] + 1) * (n1['trust'] / n1['total'])

# Rule 3: Users with more distrust than trust ratings get 1.5x negative weight
# This rule applies to N2 (distrust=3 > trust=1).
# The "1.5x negative weight" is interpreted as modifying the trust/distrust ratio calculation.
n2_modified_denominator = n2['trust'] + 1.5 * n2['distrust']
n2_contrib = -1 / (n2['total'] + 1) * (n2['trust'] / n2_modified_denominator)

# Sum all contributions to get the final score
total_score = p1_contrib + p2_contrib + p3_contrib + n1_contrib + n2_contrib

# --- Output the results step-by-step ---

print("B's Importance Score (FAT) Calculation:")
print("-" * 40)
print("Contributions from Positive Edges (Trust Incoming):")
print(f"P1: 1/({p1['total']}+1) = {p1_contrib:.4f}")
print(f"P2: 1/({p2['total']}+1) = {p2_contrib:.4f}")
print(f"P3: 1/({p3['total']}+1) = {p3_contrib:.4f}")
print("-" * 40)

print("Contributions from Negative Edges (Distrust Incoming):")
print(f"N1: -1/({n1['total']}+1) * ({n1['trust']}/{n1['total']}) = {n1_contrib:.4f}")
print(f"N2 (with 1.5x weight): -1/({n2['total']}+1) * ({n2['trust']}/({n2['trust']}+1.5*{n2['distrust']})) = {n2_contrib:.4f}")
print("-" * 40)

print("Final Equation and Score:")
print(f"Score = (1/({p1['total']}+1)) + (1/({p2['total']}+1)) + (1/({p3['total']}+1)) + (-1/({n1['total']}+1) * {n1['trust']}/{n1['total']}) + (-1/({n2['total']}+1) * {n2['trust']}/({n2['trust']}+1.5*{n2['distrust']}))")
print(f"Score = {p1_contrib:.4f} + {p2_contrib:.4f} + {p3_contrib:.4f} + {n1_contrib:.4f} + {n2_contrib:.4f}")
print(f"Final Score = {total_score:.4f}")