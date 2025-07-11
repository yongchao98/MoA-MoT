import fractions

# This script calculates the numerical coefficient of the scalar curvature term
# in the a_2 Seeley-DeWitt coefficient for a massless gauged Dirac spinor field.

# The coefficient is a sum of two contributions.

# 1. The universal contribution from the heat kernel expansion for any
#    second-order operator of the form (Laplacian - Potential).
#    This coefficient is 1/6.
c1_str = "1/6"
c1 = fractions.Fraction(c1_str)

# 2. The contribution from the potential term E in the Lichnerowicz formula
#    D^2 = Delta - E. The relevant part of E is (1/4)R.
#    This coefficient is 1/4.
c2_str = "1/4"
c2 = fractions.Fraction(c2_str)

# The total coefficient is the sum of these two parts.
total_coeff = c1 + c2

print("The numerical coefficient 'c' for the scalar curvature term in the a_2 integrand is calculated as the sum of a universal part and a part from the Lichnerowicz formula.")
print("\nThe final coefficient is the sum:")
print(f"c = {c1_str} + {c2_str}")
print(f"c = {total_coeff.numerator}/{total_coeff.denominator}")
print(f"In decimal form, c = {float(total_coeff):.4f}")
