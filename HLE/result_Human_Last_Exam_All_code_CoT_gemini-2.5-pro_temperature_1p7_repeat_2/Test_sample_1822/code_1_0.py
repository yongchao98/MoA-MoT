import sys
# Define the parameters from the problem statement.

# c: Confidence of the discriminatory pattern in the subgraph G.
c = 0.85
# s: Support of the discriminatory pattern in the subgraph G.
s = 0.12
# c_prime: Confidence of the anti-correlated pattern in the 2-hop neighborhood.
c_prime = 0.78
# beta: The bottleneck coefficient of the knowledge graph.
beta = 0.23

# The formula to calculate the minimum evidence ratio (E_r) is derived by modeling
# the problem as a signal-to-noise ratio.
# Signal (Evidence for violation) = c * s
# Noise (Conflicting evidence) = c_prime
# The bottleneck (beta) discounts the effect of the noise from the neighborhood,
# so the effective noise is c_prime * (1 - beta).
# E_r = Signal / Effective_Noise

# Calculate the numerator (Signal)
numerator = c * s
# Calculate the denominator (Effective Noise)
denominator = c_prime * (1 - beta)

# Calculate the final evidence ratio E_r
# We check for a non-zero denominator to avoid division by zero errors.
if denominator == 0:
    print("Error: Calculation resulted in division by zero.", file=sys.stderr)
    e_r = float('inf')
else:
    e_r = numerator / denominator

# Print the final equation with all the numerical values substituted in,
# as requested in the problem description.
print("Calculation of the Minimum Evidence Ratio (E_r):")
print(f"E_r = (c * s) / (c' * (1 - β))")
print(f"E_r = ({c} * {s}) / ({c_prime} * (1 - {beta}))")
print(f"E_r = {numerator} / ({c_prime} * {1-beta})")
print(f"E_r = {numerator} / {denominator}")
print(f"E_r = {e_r}")
print("\nTherefore, the minimum ratio of evidence required is approximately {:.4f}.".format(e_r))
print(f"<<<{e_r}>>>")