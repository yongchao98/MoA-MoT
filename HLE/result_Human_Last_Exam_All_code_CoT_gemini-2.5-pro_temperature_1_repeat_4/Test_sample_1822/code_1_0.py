# Define the given parameters from the problem description
c = 0.85  # Confidence of the discriminatory pattern in subgraph G
s = 0.12  # Support of the discriminatory pattern in subgraph G
c_prime = 0.78  # Confidence of the anti-correlated pattern in the 2-hop neighborhood
beta = 0.23  # Bottleneck coefficient of the knowledge graph

# Calculate the strength of the evidence for a true violation.
# This is modeled as the product of the local pattern's confidence and support.
evidence_for_violation = c * s

# Calculate the strength of the evidence for a statistical artifact.
# This is modeled as the product of the counter-evidence confidence and the bottleneck coefficient.
evidence_for_artifact = c_prime * beta

# Calculate the minimum ratio of evidence (E_r) required to confirm the violation.
# This is the ratio of the evidence for the violation to the evidence for the artifact.
E_r = evidence_for_violation / evidence_for_artifact

# Print the final equation with the numerical values and the result.
print(f"E_r = (c * s) / (c' * β)")
print(f"E_r = ({c} * {s}) / ({c_prime} * {beta})")
print(f"E_r = {evidence_for_violation} / {evidence_for_artifact}")
print(f"E_r = {E_r}")

# The final answer is the calculated value of E_r.
# Format the answer for direct extraction.
# We round the final answer to 4 decimal places for clarity.
final_answer = round(E_r, 4)
print(f"<<<{final_answer}>>>")