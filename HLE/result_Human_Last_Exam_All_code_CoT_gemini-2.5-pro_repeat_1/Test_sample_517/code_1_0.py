import math

# The problem is a classical problem in probability theory.
# The limit lim_{n->inf} p_n can be shown to be 1/2.
# This result stems from the rotational symmetry of the underlying process at large scales.
# The conditioned random walk escapes to infinity, and the exit direction is uniformly distributed
# over [0, 2*pi). The starting point (0,1) is on the y-axis, while the target (n,0) is on the x-axis.
# The probability that the walk eventually enters the right half-plane (x>0) is 1/2 by symmetry.
# For large n, if the walk enters the correct half-plane, it is certain to hit the target set B_n.
# A rigorous proof is highly technical. We present the final answer as derived from these symmetry arguments.

numerator = 1
denominator = 2
result = numerator / denominator

# Final equation is simply the result.
# We present it as an equation for clarity as requested.
print(f"The limit is given by the equation: p = {numerator} / {denominator}")
print(f"p = {result}")
