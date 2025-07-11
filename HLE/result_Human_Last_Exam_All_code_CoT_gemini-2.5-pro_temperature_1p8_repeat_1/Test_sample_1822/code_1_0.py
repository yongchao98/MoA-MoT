# Given parameters from the automated knowledge graph FAIRness auditing system
# Confidence of the discriminatory pattern P in subgraph G
c = 0.85
# Support of the discriminatory pattern P in subgraph G
s = 0.12
# Confidence of the anti-correlated pattern in the 2-hop neighborhood
c_prime = 0.78
# Bottleneck coefficient of the knowledge graph
beta = 0.23

# Step 1: Calculate the strength of the evidence for the violation (E_for).
# This is modeled as the product of confidence and support.
E_for = c * s

# Step 2: Calculate the strength of the evidence for it being a statistical artifact (E_against).
# This is modeled as the counter-evidence confidence (c_prime) amplified by the structural
# explanation (inversely proportional to beta).
E_against = c_prime / beta

# Step 3: Calculate the final evidence ratio (E_r).
# This is the ratio of E_for to E_against.
# E_r = (c * s) / (c_prime / beta) which simplifies to (c * s * beta) / c_prime
E_r = (c * s * beta) / c_prime

# Output the final result, showing the formula with the original values.
print(f"The minimum ratio of evidence (E_r) required is calculated as follows:")
print(f"E_r = (c * s * β) / c'")
print(f"E_r = ({c} * {s} * {beta}) / {c_prime}")
print(f"E_r = {E_r}")