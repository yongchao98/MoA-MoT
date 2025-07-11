import numpy as np

# These are the transition rates given in the problem
lambda01 = 0.019
lambda10 = 0.65
lambda12 = 0.4
lambda21 = 0.392
lambda23 = 0.008
lambda31 = 0.008

# At steady state (t -> +inf), the derivatives P_i'(t) are 0.
# This gives a system of linear equations. Let p_i be the steady-state probabilities P_i(+inf).
#
# 0 = -lambda01*p0 + lambda10*p1  => p0 = (lambda10 / lambda01) * p1
# 0 = lambda12*p1 - (lambda21 + lambda23)*p2 => p2 = (lambda12 / (lambda21 + lambda23)) * p1
# 0 = lambda23*p2 - lambda31*p3 => p3 = (lambda23 / lambda31) * p2
# And the normalization condition: p0 + p1 + p2 + p3 = 1
#
# Let's express p0, p2, and p3 as multiples of p1.
# Let p0 = c0*p1, p2 = c2*p1, p3 = c3*p2 = c3*c2*p1

c0 = lambda10 / lambda01
c2 = lambda12 / (lambda21 + lambda23)
c3 = lambda23 / lambda31

# Now substitute these into the normalization equation:
# c0*p1 + p1 + c2*p1 + c3*c2*p1 = 1
# p1 * (c0 + 1 + c2 + c3*c2) = 1
# p1 = 1 / (c0 + 1 + c2 + c3*c2)

denominator = c0 + 1 + c2 + c3 * c2

# Now we can calculate the values for p1 and then p0
p1 = 1 / denominator
p0 = c0 * p1

# The value we need to find is p0 + p1
result = p0 + p1

# As requested, we print the components of the final sum.
print(f"The steady-state probability P0(+inf) is: {p0}")
print(f"The steady-state probability P1(+inf) is: {p1}")
print(f"The final sum P0(+inf) + P1(+inf) is: {p0} + {p1} = {result}")