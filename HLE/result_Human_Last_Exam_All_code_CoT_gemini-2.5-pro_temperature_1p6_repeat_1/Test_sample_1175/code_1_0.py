import math
from fractions import Fraction

# Based on the analytical derivation, the expression simplifies to the form:
# (3/2) * 10^(10/3) + 37/4

# Define the components of the final equation using fractions for precision
coeff1_frac = Fraction(3, 2)
term2_frac = Fraction(37, 4)
base = 10
exponent_num = 10
exponent_den = 3

# Calculate the numerical value of the expression
result = float(coeff1_frac) * (base ** (exponent_num / exponent_den)) + float(term2_frac)

# Print the final equation with each of its numerical components, and the final calculated value.
print(f"The final simplified expression is: ({coeff1_frac.numerator}/{coeff1_frac.denominator}) * {base}^({exponent_num}/{exponent_den}) + ({term2_frac.numerator}/{term2_frac.denominator})")
print(f"The calculated value is: {result}")