import math

# Given parameters
# c: Confidence of the discriminatory pattern in subgraph G
c = 0.85
# s: Support of the discriminatory pattern in subgraph G
s = 0.12
# c_prime: Confidence of the anti-correlation in the 2-hop neighborhood
c_prime = 0.78
# beta: Bottleneck coefficient of the knowledge graph
beta = 0.23

# --- Calculation ---
# The formula for the evidence ratio (E_r) is derived by comparing the evidence for
# a true violation against the evidence for a statistical artifact.
#
# E_r = (Evidence for Violation) / (Evidence for Artifact)
# Evidence for Violation is modeled as c / (1 - s)
# Evidence for Artifact is modeled as c' / beta
#
# E_r = (c / (1 - s)) / (c_prime / beta)
# E_r = (c * beta) / (c_prime * (1 - s))

# Numerator of the final formula
numerator = c * beta
# Denominator of the final formula
denominator = c_prime * (1 - s)

# Calculate the final evidence ratio
e_r = numerator / denominator

# --- Output ---
# Print the final equation with the numbers substituted in, as requested.
print("The formula for the evidence ratio is:")
print("E_r = (c * β) / (c' * (1 - s))")
print("\nSubstituting the given values:")
print(f"E_r = ({c} * {beta}) / ({c_prime} * (1 - {s}))")
print(f"E_r = {numerator} / ({c_prime} * {1-s})")
print(f"E_r = {numerator} / {denominator}")
print("\nResult:")
print(f"The minimum ratio of evidence (E_r) required is: {e_r}")
