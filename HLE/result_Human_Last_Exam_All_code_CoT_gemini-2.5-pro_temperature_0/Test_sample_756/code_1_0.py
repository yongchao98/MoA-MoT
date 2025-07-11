import numpy as np

# Let the quadratic be f(x) = ax^2 + bx + c.
# We are given |f(x)| <= 1 for all x in [-1, 1].
# We want to find the maximum value of |b| + |c|.

# Based on the analysis, the maximum value is achieved for polynomials
# where the vertex lies on the boundary of the interval [-1, 1].
# One such polynomial is f(x) = -0.5*x^2 + 1*x + 0.5.
# Let's verify this polynomial.
# The vertex is at x = -1/(2*(-0.5)) = 1.
# f(1) = -0.5 + 1 + 0.5 = 1.
# At the other end of the interval, f(-1) = -0.5 - 1 + 0.5 = -1.
# Since the parabola opens downwards and its vertex is at x=1,
# for x in [-1, 1], the function is monotonic.
# The values f(x) are between f(-1)=-1 and f(1)=1.
# So, |f(x)| <= 1 is satisfied.

# For this polynomial, the coefficients are:
a = -0.5
b = 1
c = 0.5

# The value of |b| + |c| is:
max_value = abs(b) + abs(c)

print(f"An optimal polynomial is f(x) = ({a})x^2 + ({b})x + ({c}).")
print("For this polynomial, the value of |b| + |c| is:")
print(f"|{b}| + |{c}| = {max_value}")
