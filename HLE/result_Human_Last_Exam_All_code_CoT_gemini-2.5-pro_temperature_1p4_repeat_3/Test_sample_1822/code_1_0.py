import math

# Define the given parameters for the fairness audit
# c: Confidence of the discriminatory pattern P in subgraph G
c = 0.85
# s: Support for the pattern P. This is contextual information about the pattern's
# prevalence and is not used in the evidence ratio calculation itself.
s = 0.12
# c_prime: Confidence of the anti-correlated pattern in the neighborhood
c_prime = 0.78
# beta: The bottleneck coefficient of the knowledge graph
beta = 0.23

# The formula for the minimum evidence ratio (E_r) is derived by modeling the
# ratio of effective confirmatory evidence to confounding evidence.
#
# Effective Confirmatory Evidence = c * beta
# The initial confidence 'c' is discounted by the bottleneck coefficient 'beta'.
# A low beta indicates a strong bottleneck, increasing the chance of statistical
# artifacts and thus reducing the effective strength of the evidence.
#
# Confounding Evidence = c_prime
# The confidence of the contradictory finding is the direct confounding evidence.
#
# E_r = (c * beta) / c_prime

# Calculate the numerator (effective confirmatory evidence)
numerator = c * beta

# Calculate the final evidence ratio E_r
E_r = numerator / c_prime

# Print the final equation with all the numbers, as requested.
print("Calculating the Evidence Ratio (E_r):")
print(f"E_r = (c * beta) / c_prime")
print(f"E_r = ({c} * {beta}) / {c_prime}")
print(f"E_r = {numerator} / {c_prime}")
print(f"E_r = {E_r}")
