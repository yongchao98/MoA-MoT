# Define the given parameters
c = 0.85  # Confidence of the discriminatory pattern P
s = 0.12  # Support of the discriminatory pattern P
c_prime = 0.78  # Confidence of the anti-correlation in the 2-hop neighborhood
beta = 0.23  # Bottleneck coefficient of the knowledge graph

# Step 1: Calculate the evidence for the discriminatory pattern.
# This is the product of its confidence and support.
evidence_for = c * s

# Step 2: Calculate the effective counter-evidence from the neighborhood.
# This is the confidence of the anti-correlation, modulated by the graph's structure.
# (1 - beta) represents the influence of the neighborhood.
effective_counter_evidence = c_prime * (1 - beta)

# Step 3: Calculate the minimum ratio of evidence (E_r).
# This is the ratio of the supporting evidence to the counter-evidence.
E_r = evidence_for / effective_counter_evidence

# Print the final equation and the result
print(f"The formula for the minimum ratio of evidence (E_r) is:")
print(f"E_r = (c * s) / (c' * (1 - β))")
print(f"\nPlugging in the values:")
print(f"E_r = ({c} * {s}) / ({c_prime} * (1 - {beta}))")
print(f"E_r = {evidence_for} / ({c_prime} * {1 - beta})")
print(f"E_r = {evidence_for} / {effective_counter_evidence}")
print(f"\nThe minimum ratio of evidence required is: {E_r}")
