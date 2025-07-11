# Define the given parameters from the problem description
c = 0.85  # Confidence of the initial discriminatory pattern in subgraph G
s = 0.12  # Support of the discriminatory pattern in subgraph G
c_prime = 0.78  # Confidence of the anti-correlation in the 2-hop neighborhood
beta = 0.23  # Bottleneck coefficient of the knowledge graph

# The formula to calculate the minimum ratio of evidence (E_r) is derived by
# comparing the structurally adjusted evidence for the violation against the
# evidence for it being a statistical artifact.
#
# E_r = (Strength of Evidence FOR Violation) / (Strength of Evidence AGAINST Violation)
#
# The evidence for the violation is c * s.
# The bottleneck coefficient (beta) acts as a discounting factor on local evidence. A smaller beta
# signifies a stronger bottleneck, making local patterns more likely to be artifacts.
# So, the structurally adjusted evidence for the violation is c * s * beta.
# The evidence against the violation is the confidence of the counter-evidence, c_prime.
#
# Final Formula: E_r = (c * s * beta) / c_prime

# Calculate the numerator (structurally adjusted evidence for violation)
numerator = c * s * beta

# The denominator is the evidence against the violation
denominator = c_prime

# Calculate the final evidence ratio E_r
E_r = numerator / denominator

# Print the breakdown of the calculation and the final result
print("To determine the minimum ratio of evidence (E_r) required to confirm a true violation, we use the following formula:")
print("E_r = (c * s * β) / c'")
print("\nPlugging in the given values:")
print(f"E_r = ({c} * {s} * {beta}) / {c_prime}")
print(f"E_r = {numerator:.5f} / {denominator}")
print(f"\nThe calculated minimum ratio of evidence is:")
print(f"E_r = {E_r}")

print(f"\n<<<{E_r}>>>")