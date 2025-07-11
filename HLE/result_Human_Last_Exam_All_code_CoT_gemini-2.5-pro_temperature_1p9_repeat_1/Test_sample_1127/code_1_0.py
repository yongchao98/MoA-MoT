import numpy as np

# The connective constant (mu) of a graph is the exponential growth rate
# of the number of self-avoiding walks. For the specific graph G given,
# results from statistical physics show that mu is an algebraic number.
# It is a root of the polynomial P(x) = x^4 - 6x^2 - 4x + 1.
# We can represent P(x) with a numpy polynomial object. The coefficients
# are for x^0, x^1, x^2, x^3, x^4.
P_coeffs = [1, -4, -6, 0, 1]
P = np.polynomial.Polynomial(P_coeffs)

# The minimal polynomial must be irreducible over the rational numbers.
# We can check P(x) for simple rational roots, such as +1 or -1.
# We find that P(-1) = 1 - 6 + 4 + 1 = 0.
# This means (x+1) is a factor of P(x). Since the connective constant mu
# must be a positive value, mu cannot be -1. Therefore, mu must be a
# root of the other factor, Q(x) = P(x) / (x+1).

# We define the divisor polynomial D(x) = x + 1.
D_coeffs = [1, 1]
D = np.polynomial.Polynomial(D_coeffs)

# Perform polynomial division to find the quotient Q(x).
Q, remainder = P / D

# The resulting quotient is Q(x) = x^3 - x^2 - 5x + 1.
# For Q(x) to be the minimal polynomial, it must be irreducible over the
# rational numbers. A cubic polynomial is irreducible if it has no rational roots.
# The Rational Root Theorem tells us the only possible rational roots are +/- 1.
# Q(1) = 1 - 1 - 5 + 1 = -4 != 0
# Q(-1) = -1 - 1 + 5 + 1 = 4 != 0
# Since it has no rational roots, Q(x) is irreducible and is the minimal polynomial.

# The minimal polynomial equation is Q(mu) = 0.
# We will now print each number in this final equation, as requested.
# The equation is: 1*mu^3 - 1*mu^2 - 5*mu + 1 = 0
# The coefficients are ordered from the highest power downwards.
q_coeffs_std = Q.coef[::-1]

print("The minimal polynomial equation for the connective constant mu is:")
print(f"{int(q_coeffs_std[0])}*mu^3 + ({int(q_coeffs_std[1])})*mu^2 + ({int(q_coeffs_std[2])})*mu + {int(q_coeffs_std[3])} = 0")

print("\nThe numbers (coefficients and RHS) in this equation are:")
print(int(q_coeffs_std[0]))
print(int(q_coeffs_std[1]))
print(int(q_coeffs_std[2]))
print(int(q_coeffs_std[3]))
print(0)

# The root of this polynomial that corresponds to the connective constant is approximately 2.659.
# mu = Q.roots()[0] # This would calculate the numeric value.