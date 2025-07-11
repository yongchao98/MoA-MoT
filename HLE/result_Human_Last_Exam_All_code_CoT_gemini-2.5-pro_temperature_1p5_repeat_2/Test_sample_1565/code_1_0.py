import sys

# Define a function to handle potential division by zero
def safe_div(numerator, denominator):
    if denominator == 0:
        return 0
    return numerator / denominator

# --- User Data ---
# Trust Incoming
# P1: 7 total, 4 trust-senior, 2 trust-peer, 1 distrust-peer
p1_total = 7
p1_trust = 4 + 2
p1_distrust = 1

# P2: 6 total, 3 trust-junior, 1 trust-senior, 2 distrust-peer
p2_total = 6
p2_trust = 3 + 1
p2_distrust = 2

# P3: 4 total, 2 trust-peer, 1 distrust-senior, 1 distrust-junior
p3_total = 4
p3_trust = 2
p3_distrust = 1 + 1

# Distrust Incoming
# N1: 6 total, 3 trust-senior, 2 distrust-peer, 1 distrust-junior
n1_total = 6
n1_trust = 3
n1_distrust = 2 + 1

# N2: 4 total, 1 trust-peer, 3 distrust-junior
n2_total = 4
n2_trust = 1
n2_distrust = 3

# --- Calculations ---

# Rule 1: Positive edge contributes 1/(total_relationships + 1)
p1_score = 1 / (p1_total + 1)
p2_score = 1 / (p2_total + 1)
p3_score = 1 / (p3_total + 1)

# Rule 2: Negative edge with mixed ratings: -1 /(total relationships+1) x (trust_ratings/total)
# For N1
n1_score = -1 / (n1_total + 1) * safe_div(n1_trust, n1_total)
# Rule 3 Check for N1: Users with more distrust than trust ratings get 1.5x negative weight
if n1_distrust > n1_trust:
    n1_score *= 1.5

# For N2
n2_score = -1 / (n2_total + 1) * safe_div(n2_trust, n2_total)
# Rule 3 Check for N2
if n2_distrust > n2_trust:
    n2_score *= 1.5

# Sum the scores for the final FAT measure
total_score = p1_score + p2_score + p3_score + n1_score + n2_score

# --- Output ---
# Redirect print to stdout for web UI
sys.stdout.reconfigure(encoding='utf-8')

print("Calculating B's FAT score:")
print("Contribution from P1 (Trust): 1 / (7 + 1)")
print("Contribution from P2 (Trust): 1 / (6 + 1)")
print("Contribution from P3 (Trust): 1 / (4 + 1)")
print("Contribution from N1 (Distrust): -1 / (6 + 1) * (3 / 6)")
print("Contribution from N2 (Distrust): [-1 / (4 + 1) * (1 / 4)] * 1.5 (since distrust 3 > trust 1)")
print("\nFinal Equation:")
print(f"{p1_score:.4f} (from P1) + {p2_score:.4f} (from P2) + {p3_score:.4f} (from P3) + ({n1_score:.4f}) (from N1) + ({n2_score:.4f}) (from N2) = {total_score:.4f}")
