import numpy as np
from scipy.special import zeta
from math import factorial
from fractions import Fraction

# Plan: Calculate the final result by summing the four integral components
# as derived in the step-by-step analysis.

# I_1 = pi^2/6
c1_num_pi2 = 1
c1_den_pi2 = 6

# I_2 = (8/15)*pi^8
# This comes from Gamma(8)*zeta(8) = 7! * pi^8/9450
c2_frac = Fraction(factorial(7), 9450)
c2_num_pi8 = c2_frac.numerator
c2_den_pi8 = c2_frac.denominator

# I_3 = pi^2/6 - 1
c3_num_pi2 = 1
c3_den_pi2 = 6
c3_const = -1

# I_4 = 2 - pi/2
c4_const = 2
c4_num_pi = -1
c4_den_pi = 2

# Summing the coefficients for each term (pi^8, pi^2, pi, constant)
# Coefficient for pi^8 term:
c_pi8 = Fraction(c2_num_pi8, c2_den_pi8)

# Coefficient for pi^2 term:
c_pi2 = Fraction(c1_num_pi2, c1_den_pi2) + Fraction(c3_num_pi2, c3_den_pi2)

# Coefficient for pi term:
c_pi = Fraction(c4_num_pi, c4_den_pi)

# Coefficient for constant term:
c_const = c3_const + c4_const

print("The value of the integral is of the form: C1 * pi^8 + C2 * pi^2 + C3 * pi + C4")
print("The coefficients are calculated as follows:")
print(f"C1 (for pi^8) = {c_pi8.numerator}/{c_pi8.denominator}")
print(f"C2 (for pi^2) = {c_pi2.numerator}/{c_pi2.denominator}")
print(f"C3 (for pi) = {c_pi.numerator}/{c_pi.denominator}")
print(f"C4 (constant) = {c_const}")

print("\nThe final equation is:")
print(f"{c_pi8.numerator}/{c_pi8.denominator} * pi^8 + {c_pi2.numerator}/{c_pi2.denominator} * pi^2 - {abs(c_pi.numerator)}/{abs(c_pi.denominator)} * pi + {c_const}")
