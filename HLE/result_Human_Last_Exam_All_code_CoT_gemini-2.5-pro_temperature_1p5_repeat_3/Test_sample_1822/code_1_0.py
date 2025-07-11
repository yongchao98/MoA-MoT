# Define the given parameters from the knowledge graph FAIRness auditing system
# Confidence of the discriminatory pattern P in subgraph G
c = 0.85
# Support of the pattern P in subgraph G
s = 0.12
# Confidence of the anti-correlated pattern in the 2-hop neighborhood of G
c_prime = 0.78
# Bottleneck coefficient based on the hourglass structural principle
beta = 0.23

# --- Calculation of the Evidence Ratio (E_r) ---

# The evidence for the violation is based on the confidence (c) and support (s)
# of the pattern in the local subgraph.
evidence_for_base = c * s

# This local evidence is discounted by the bottleneck coefficient (beta),
# as a strong bottleneck makes local patterns less reliable.
adjusted_evidence_for = evidence_for_base * beta

# The evidence against the violation is the confidence of the contradictory
# pattern in the wider, more stable 2-hop neighborhood.
evidence_against = c_prime

# The final ratio E_r is the adjusted evidence for the violation
# divided by the evidence against it.
E_r = adjusted_evidence_for / evidence_against

# --- Output the results ---

# Print the formula with the substituted values
print(f"The formula for the evidence ratio E_r is (c * s * beta) / c_prime")
print(f"E_r = ({c} * {s} * {beta}) / {c_prime}")

# Print the final calculated value
print(f"The minimum ratio of evidence (E_r) required is: {E_r}")
<<<0.03007692307692308>>>