import math

# The problem is to find the asymptotic value of d_{B,delta}, which represents
# the minimum degree of a polynomial p satisfying |p(x) - e^{-x}| < delta
# for all x in [0, B].
# We use the notation L = log(delta^{-1}).

# The analysis, based on approximation theory, shows two distinct asymptotic behaviors
# depending on the ratio L/B:
# 1. If L/B is small, d is proportional to sqrt(B*L).
# 2. If L/B is large, d is proportional to L / log(L/B).

# A single formula that combines these behaviors and provides the asymptotic
# value for d_{B,delta} is constructed. The code below prints this formula.
# In the formula, 'sqrt' is the square root function and 'log' is the natural logarithm.
# The constant 2.718281828 is an approximation of Euler's number 'e',
# used to ensure the argument of the logarithm is always greater than 1.

# The final asymptotic expression for d is:
final_equation = "d = sqrt(B*L) + L / log(L/B + 2.718281828)"

# The numbers in this equation are 1 (implicitly for B, L), and e.
print(final_equation)
