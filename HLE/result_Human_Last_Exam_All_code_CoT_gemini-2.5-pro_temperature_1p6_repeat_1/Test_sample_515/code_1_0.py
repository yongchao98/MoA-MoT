import math

# The problem asks for the probability that a 2D simple random walk,
# starting at x0 = (0,1) and conditioned to never enter the origin,
# eventually reaches the set of neighbors of y = (3600,0).

# The derivation in the thought process shows that this probability p can be approximated by:
# p = 1/2 + (a(y) - a(y-x0)) / (2 * a(x0))
# where a(x) is the potential kernel of the 2D SRW.

# Using x0 = (0,1) and y = (3600,0)
# a(x0) = a(0,1) = 1/pi
# The difference a(y) - a(y-x0) is approximated by -1/(pi * 3600^2)

# So, p approx 1/2 + (-1/(pi * 3600^2)) / (2 * (1/pi))
# p approx 1/2 - 1 / (2 * 3600^2)

Y = 3600
# The equation for the probability p is derived as:
# p = 1/2 - 1/(2*Y^2)
prob = 0.5 - 1 / (2 * Y**2)

print(f"The starting point is (0,1). The target point y is ({Y},0).")
print(f"The equation for the probability is approximately p = 1/2 - 1/(2 * Y^2)")
print(f"p = 1/2 - 1/(2 * {Y}^2)")
print(f"p = {0.5} - 1/({2 * Y**2})")
print(f"p = {prob}")
print(f"The approximate answer with two significant digits is {prob:.2f}")

# Return the final numerical answer rounded to two significant digits.
final_answer = round(prob, 2)