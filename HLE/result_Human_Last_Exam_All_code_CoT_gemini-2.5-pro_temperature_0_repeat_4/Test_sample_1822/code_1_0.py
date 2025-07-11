# Define the given parameters
# c: Confidence of the discriminatory pattern in the subgraph G
c = 0.85
# s: Support of the discriminatory pattern in the subgraph G
s = 0.12
# c_prime: Confidence of the anti-correlation pattern in the 2-hop neighborhood
c_prime = 0.78
# beta: Bottleneck coefficient of the knowledge graph
beta = 0.23

# The formula for the ratio of evidence (E_r) is derived by comparing
# the evidence for a true violation against the evidence for a statistical artifact.
# Evidence for Violation = c * s
# Evidence for Artifact = c_prime / beta
# E_r = (Evidence for Violation) / (Evidence for Artifact) = (c * s * beta) / c_prime

# Calculate the numerator (Evidence for Violation, adjusted by bottleneck factor)
numerator = c * s * beta

# Calculate the final ratio E_r
E_r = numerator / c_prime

# Print the equation with the numbers plugged in, as requested.
print("Step 1: Define the formula for the evidence ratio (E_r)")
print("E_r = (c * s * β) / c'")
print("\nStep 2: Substitute the given values into the formula")
print(f"E_r = ({c} * {s} * {beta}) / {c_prime}")
print(f"E_r = ({c * s} * {beta}) / {c_prime}")
print(f"E_r = {numerator} / {c_prime}")
print("\nStep 3: Calculate the final result")
print(f"E_r = {E_r}")
