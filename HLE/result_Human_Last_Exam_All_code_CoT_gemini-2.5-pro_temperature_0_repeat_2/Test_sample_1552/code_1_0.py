from fractions import Fraction

# The problem is to find the numerical coefficient C in the integrand of the second heat kernel coefficient a_2.
# This coefficient appears in the term C * R * tr(I), where R is the scalar curvature and I is the identity.
# The formula for C is derived from two components:
# 1. A universal geometric term, which contributes a factor of 1/6.
# 2. A term from the Lichnerowicz formula for the squared Dirac operator (D^2), which contributes a factor of 1/4.
# The final coefficient is the difference between these two values, as dictated by the heat kernel expansion formula.

universal_term_coeff = Fraction(1, 6)
lichnerowicz_term_coeff = Fraction(1, 4)

# The formula for the coefficient is (universal_term - lichnerowicz_term)
final_coeff = universal_term_coeff - lichnerowicz_term_coeff

print("The final coefficient is calculated by subtracting the Lichnerowicz formula contribution from the universal geometric term.")
print("The equation for the coefficient is:")
print(f"{universal_term_coeff.numerator}/{universal_term_coeff.denominator} - {lichnerowicz_term_coeff.numerator}/{lichnerowicz_term_coeff.denominator} = {final_coeff.numerator}/{final_coeff.denominator}")