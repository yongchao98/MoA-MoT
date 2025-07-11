import math
from fractions import Fraction

# The problem is solved symbolically. This code will format and print the final derived expression.
# The step-by-step derivation is outlined above.

# The final ratio is composed of coefficients for pi, sqrt(3), and a rational term.
pi_coeff = Fraction(1, 6)
sqrt3_coeff = Fraction(1, 4)
rational_coeff = Fraction(-1, 4)

# Print the final result in its exact form, showing each number in the equation.
print("The ratio of the area of D to the area of S is:")
print(f"pi/{pi_coeff.denominator} + sqrt(3)/{sqrt3_coeff.denominator} - {abs(rational_coeff.numerator)}/{abs(rational_coeff.denominator)}")