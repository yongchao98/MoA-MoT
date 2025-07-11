import math

# Given parameters
epsilon = 0.05  # Marginal completeness tolerance
confidence_level = 0.99
alpha = 2.5  # Shape parameter for truncated Pareto distribution
gamma = 2.1  # Power-law exponent for scale-free properties

# The problem asks for a minimum sampling ratio 'r'.
# Standard sample size formulas calculate an absolute number 'n', and the ratio r = n/N
# would depend on the total number of triples N, which is unknown.
# The presence of specific graph complexity parameters (alpha, gamma) suggests
# a model where the required ratio depends on the graph's structural complexity.
# A heuristic model relates the ratio 'r' to the tolerance 'epsilon' and a
# complexity factor derived from the exponents alpha and gamma.
# The terms (gamma - 1) and (alpha - 1) are related to the moments of these distributions.

# Formula for the sampling ratio r
# r = epsilon / ((gamma - 1) + (alpha - 1))
r_numerator = epsilon
r_denominator = (gamma - 1) + (alpha - 1)
r = r_numerator / r_denominator

# Round the final result to 4 decimal places
r_rounded = round(r, 4)

# Print the equation with the final numbers
print(f"The calculation uses the formula: r = ε / ((γ - 1) + (α - 1))")
print(f"Plugging in the values:")
print(f"r = {epsilon} / (({gamma} - 1) + ({alpha} - 1))")
print(f"r = {epsilon} / ({gamma - 1} + {alpha - 1})")
print(f"r = {epsilon} / {r_denominator}")
print(f"r ≈ {r:.6f}")
print(f"The minimum ratio r, rounded to 4 decimal places, is: {r_rounded}")

# Final Answer
# print(f'<<<r={r_rounded}>>>')