# Define the given parameters
# c: Confidence of the discriminatory pattern in the subgraph G
c = 0.85
# s: Support of the discriminatory pattern in the subgraph G
s = 0.12
# c_prime: Confidence of the anti-correlated pattern in the 2-hop neighborhood
c_prime = 0.78
# beta: Bottleneck coefficient of the knowledge graph
beta = 0.23

# The problem asks for the minimum ratio of evidence (E_r) required to confirm
# the violation is true and not a statistical artifact.

# We model the evidence for the violation as its confidence 'c', but discounted by the
# bottleneck coefficient 'beta'. A strong bottleneck (low beta) makes local
# observations less reliable.
# Adjusted Evidence for Violation = c * beta

# The evidence against the violation (i.e., for it being an artifact) is the
# confidence of the contradictory finding 'c_prime'.

# The ratio E_r is therefore the adjusted evidence for divided by the evidence against.
E_r = (c * beta) / c_prime

# Print the final equation with all the numbers and the result.
print(f"The minimum ratio of evidence (E_r) is calculated as:")
print(f"E_r = (c * beta) / c_prime")
print(f"E_r = ({c} * {beta}) / {c_prime}")
print(f"E_r = {c*beta} / {c_prime}")
print(f"E_r = {E_r}")
