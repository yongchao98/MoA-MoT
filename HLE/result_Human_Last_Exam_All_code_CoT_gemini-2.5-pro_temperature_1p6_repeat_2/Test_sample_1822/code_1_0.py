import sys

# Given parameters from the problem description
# c: Confidence of the discriminatory pattern in the subgraph G
c = 0.85
# s: Support of the discriminatory pattern in the subgraph G
s = 0.12
# c_prime: Confidence of the anti-correlated pattern in the 2-hop neighborhood
c_prime = 0.78
# beta: Bottleneck coefficient of the knowledge graph
beta = 0.23

# --- Step 1: Model the evidence FOR the violation ---
# The evidence for the violation is the product of its local confidence (c) and support (s).
# This is scaled by the bottleneck coefficient (beta), which adjusts the significance
# of local findings based on the graph's structure. A strong bottleneck (low beta)
# reduces the weight of local evidence.
evidence_pro = c * s * beta

# --- Step 2: Model the evidence AGAINST the violation ---
# The evidence against the violation is the confidence of the contradictory pattern (c_prime)
# found in the wider neighborhood.
evidence_con = c_prime

# --- Step 3: Calculate the ratio of evidence (E_r) ---
# E_r is the ratio of the evidence for the violation to the evidence against it.
if evidence_con == 0:
    print("Error: Division by zero. The confidence of the anti-correlated pattern (c') cannot be zero.", file=sys.stderr)
    E_r = float('inf')
else:
    E_r = evidence_pro / evidence_con

# --- Step 4: Output the results ---
# The prompt requires showing the inputs and the equation used.
print("--- Input Values ---")
print(f"Confidence of discriminatory pattern (c): {c}")
print(f"Support of discriminatory pattern (s): {s}")
print(f"Confidence of anti-correlated pattern (c'): {c_prime}")
print(f"Bottleneck coefficient (β): {beta}")
print("\n--- Calculation ---")
print("The ratio of evidence (E_r) is calculated using the formula: (c * s * β) / c'")
# Here we output each number in the final equation as requested.
print(f"E_r = ({c} * {s} * {beta}) / {c_prime}")
print(f"E_r = {evidence_pro} / {evidence_con}")
print(f"\nThe minimum ratio of evidence (E_r) is: {E_r}")