# Problem parameters
d = 2  # Dimension of the ambient space R^d
alpha = 8 / 5  # Frostman exponent of the measure

# The spherically averaged Stein-Tomas inequality gives a bound on the decay exponent.
# The bound holds for an s-energy measure, where s must be in the range [0, d-1].
# Our measure is an alpha-Frostman measure, so its s-energy is finite for any s < alpha.
# We choose the largest possible s allowed by the theorem to get the tightest bound.
s = 1.0  # s must be <= d-1 = 1. This choice is valid since s < alpha (1 < 1.6).

# Calculate the exponent c from the theorem.
c = (d - 1 - s) / 2

# This calculated value represents the upper bound for the decay exponent.
# It has been shown that for dimensions alpha > 1 (like our 1.6), there exist measures
# for which this bound is sharp (i.e., the L2 norm does not decay to 0).
# This establishes that c cannot be lower than this value.
# Thus, the calculated c is the smallest possible value.

print(f"The dimension of the space is d = {d}.")
print(f"The Frostman exponent is alpha = {alpha}.")
print(f"To find the decay exponent, we use the spherically averaged Stein-Tomas inequality.")
print(f"We can choose the parameter s=1 since 1 < alpha and 1 <= d-1.")
print(f"The formula for the exponent is c = (d - 1 - s) / 2.")
print(f"Plugging in the values:")
print(f"c = ({d} - {s} - 1) / 2 = {int(c)}")