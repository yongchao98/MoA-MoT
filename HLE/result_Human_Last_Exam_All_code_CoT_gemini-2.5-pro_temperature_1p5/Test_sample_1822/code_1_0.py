# Define the given parameters
c = 0.85  # Confidence of discriminatory pattern in subgraph G
s = 0.12  # Support of discriminatory pattern in subgraph G
c_prime = 0.78  # Confidence of anti-correlated pattern in the neighborhood
beta = 0.23  # Bottleneck coefficient of the knowledge graph

# Calculate the evidence for the violation
# This is modeled as the product of the confidence and support of the pattern.
evidence_pro = c * s

# Calculate the evidence for the statistical artifact hypothesis
# This is modeled as the product of the counter-evidence confidence
# and a structural factor (1 - beta).
structural_factor = 1 - beta
evidence_con = c_prime * structural_factor

# Calculate the minimum ratio of evidence (E_r)
# E_r is the ratio of the pro-violation evidence to the con-violation evidence.
e_r = evidence_pro / evidence_con

# Print the calculation steps and the final result
print(f"The formula for the evidence ratio (E_r) is: (c * s) / (c' * (1 - β))")
print(f"E_r = ({c} * {s}) / ({c_prime} * (1 - {beta}))")
print(f"E_r = {evidence_pro} / ({c_prime} * {structural_factor})")
print(f"E_r = {evidence_pro} / {evidence_con}")
print(f"The minimum ratio of evidence (E_r) is: {e_r}")