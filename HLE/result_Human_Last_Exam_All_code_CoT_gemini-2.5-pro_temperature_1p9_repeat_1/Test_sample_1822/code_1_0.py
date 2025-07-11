# Given values from the problem description
c = 0.85  # Confidence of the discriminatory pattern in subgraph G
s = 0.12  # Support of the discriminatory pattern in subgraph G
c_prime = 0.78  # Confidence of the anti-correlated pattern in the 2-hop neighborhood
beta = 0.23  # Bottleneck coefficient of the knowledge graph

# --- Explanation of the formula ---
# The minimum ratio of evidence required (E_r) is calculated by comparing the
# weighted counter-evidence to the weighted primary evidence.
#
# E_r = (Weighted Counter-Evidence) / (Weighted Primary Evidence)
#
# - The primary evidence is the confidence (c) scaled by its support (s).
# - The counter-evidence is the confidence of the anti-correlated pattern (c').
# - The bottleneck coefficient (beta) determines the weights. A low beta means
#   the broader context is more reliable than the local pattern.
#   - Weight for local evidence (c * s) is beta.
#   - Weight for broader context evidence (c') is (1 - beta).
#
# This results in the formula: E_r = (c' * (1 - beta)) / (c * s * beta)
# A ratio > 1 suggests the counter-evidence is stronger.

# --- Calculation Step by Step ---

# 1. Calculate the weighted primary evidence (the denominator)
# This is the evidence for the violation, adjusted for its low support and
# the leaky graph structure that makes local findings less reliable.
weighted_primary_evidence = c * s * beta

# 2. Calculate the weighted counter-evidence (the numerator)
# This is the evidence against the violation, weighted by the reliability of
# the broader context.
weighted_counter_evidence = c_prime * (1 - beta)

# 3. Calculate the final evidence ratio (E_r)
if weighted_primary_evidence > 0:
    E_r = weighted_counter_evidence / weighted_primary_evidence
else:
    E_r = float('inf') # Avoid division by zero

# --- Output the results ---

print("Calculating the Minimum Ratio of Evidence Required (E_r)")
print("Formula: E_r = (c' * (1 - β)) / (c * s * β)")
print("-" * 50)
print(f"Plugging in the values:")
print(f"E_r = ({c_prime} * (1 - {beta})) / ({c} * {s} * {beta})")
print(f"E_r = ({c_prime} * {1 - beta}) / ({c * s} * {beta})")
print(f"E_r = {weighted_counter_evidence:.4f} / {weighted_primary_evidence:.4f}")
print("-" * 50)
print(f"The minimum ratio of evidence required (E_r) is: {E_r:.3f}")
print("This means the weighted counter-evidence is over 25 times stronger than the weighted primary evidence for the violation.")
