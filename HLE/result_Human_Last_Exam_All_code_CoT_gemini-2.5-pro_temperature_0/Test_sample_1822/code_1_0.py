# Define the given parameters
# c: Confidence of the discriminatory pattern in the subgraph G
c = 0.85
# s: Support of the discriminatory pattern in the subgraph G
s = 0.12
# c_prime: Confidence of the anti-correlated pattern in the 2-hop neighborhood
c_prime = 0.78
# beta: Bottleneck coefficient of the knowledge graph
beta = 0.23

# The formula for the minimum evidence ratio (E_r) is derived by comparing
# the structurally-discounted evidence for the violation against the counter-evidence.
# Evidence for violation = c * s
# The bottleneck coefficient (beta) acts as a discount factor for local evidence,
# as a strong bottleneck (low beta) makes local findings less reliable.
# Structurally-discounted evidence for violation = c * s * beta
# Counter-evidence = c_prime
# E_r = (c * s * beta) / c_prime

# Calculate the numerator (structurally-discounted evidence for violation)
numerator = c * s * beta

# Calculate the final evidence ratio
E_r = numerator / c_prime

# Output the formula and the calculation steps
print("The formula for the evidence ratio E_r is:")
print("E_r = (c * s * β) / c'")
print("\nSubstituting the given values:")
print(f"E_r = ({c} * {s} * {beta}) / {c_prime}")
print(f"E_r = {numerator} / {c_prime}")
print(f"\nThe calculated minimum ratio of evidence is:")
print(f"E_r = {E_r}")

# The final answer is the numerical value of E_r
# print(f"\n<<<{E_r}>>>")