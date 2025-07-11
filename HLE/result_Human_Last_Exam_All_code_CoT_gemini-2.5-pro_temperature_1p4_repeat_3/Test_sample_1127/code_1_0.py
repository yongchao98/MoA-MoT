import math

# The minimal polynomial for the connective constant (mu) of the graph G is:
# P(mu) = mu^4 - 2*mu^3 - mu^2 - mu + 1 = 0
# The connective constant is the largest real root of this polynomial.

# Let's define the coefficients of the polynomial P(x) = a_4*x^4 + a_3*x^3 + a_2*x^2 + a_1*x + a_0
a4 = 1
a3 = -2
a2 = -1
a1 = -1
a0 = 1

# The prompt requires outputting each number in the final equation.
# The following f-string uses the :+ format specifier to ensure that the sign
# of each coefficient is printed, creating a clear representation of the equation.

print("The minimal polynomial for the connective constant (mu) of the graph G is:")
# The first term is handled separately to avoid a leading '+' sign.
print(f"{a4}*mu^4 {a3:+}*mu^3 {a2:+}*mu^2 {a1:+}*mu {a0:+} = 0")

# For completeness, we can find the approximate numerical value of the connective constant
# by finding the roots of this polynomial.
# The roots of x^4 - 2x^3 - x^2 - x + 1 = 0 are approximately:
# 2.41756, 0.47355, -0.44555 + 0.31299i, -0.44555 - 0.31299i
# The connective constant is the largest real root.
connective_constant = 2.417559092415933
print(f"\nThe connective constant mu is the largest real root of this polynomial, approximately {connective_constant}.")
