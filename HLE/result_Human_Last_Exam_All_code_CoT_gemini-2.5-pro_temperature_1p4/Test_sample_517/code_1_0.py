import math

# This script calculates the limit of p_n by analyzing the asymptotic scaling of the probability.
# The probability p_n scales as n raised to a power, which we calculate first.

# 1. Theoretical exponents
# The transverse fluctuation exponent for a 2D random walk conditioned to avoid the origin
# (which is in the LERW/SLE2 universality class) is zeta.
zeta = 3/4

# The walk must travel a distance 'n' to reach the target region.
# The distance scales with n to the power of 1.
distance_exponent = 1

# 2. Scaling of the probability p_n
# For the walk to hit the target, its random asymptotic direction must fall within a narrow
# angular window. The size of this window is proportional to (fluctuation_width / distance).
# Fluctuation width at distance n scales as n^zeta.
# Distance is n.
# So, the angular window size scales as n^zeta / n^1 = n^(zeta - 1).
# The probability p_n is proportional to this window size.
pn_exponent = zeta - distance_exponent

# 3. Calculate the limit
# We need to compute lim_{n->inf} p_n, which is lim_{n->inf} C * n^(pn_exponent).
# Since pn_exponent is negative, this limit is 0.
final_limit = 0

# The final equation is lim_{n->inf} p_n = 0.
# The numbers involved in reaching this conclusion are the exponents.
print(f"The fluctuation exponent of the walk is zeta = {zeta}")
print(f"The probability p_n scales like n to the power of (zeta - 1).")
print(f"Calculated exponent for n: {zeta} - {distance_exponent} = {pn_exponent}")
print("\nSince the exponent is negative, the probability p_n approaches 0 as n goes to infinity.")
print("The final equation is: lim_{n->inf} p_n = 0")
print("The number in this final equation is:")
print(final_limit)
