# Define the given parameters
c = 0.85  # Confidence of the discriminatory pattern in subgraph G
s = 0.12  # Support of the discriminatory pattern in subgraph G
c_prime = 0.78  # Confidence of the anti-correlated pattern in the neighborhood
beta = 0.23  # Bottleneck coefficient of the knowledge graph

# Step 1: Calculate the supported confidence for the primary pattern (evidence for violation)
evidence_for_violation = c * s

# Step 2: Calculate the dampened confidence for the conflicting pattern (evidence against violation)
evidence_against_violation = c_prime * beta

# Step 3: Calculate the minimum ratio of evidence required
E_r = evidence_for_violation / evidence_against_violation

# Output the full equation and the final result
print(f"The minimum ratio of evidence (E_r) is calculated as:")
print(f"E_r = (c * s) / (c' * β)")
print(f"E_r = ({c} * {s}) / ({c_prime} * {beta})")
print(f"E_r = {evidence_for_violation} / {evidence_against_violation}")
print(f"E_r = {E_r}")
