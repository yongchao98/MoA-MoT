from fractions import Fraction

# Based on the analysis, the imaginary part of the sum of the integrals is:
# Im(I) = -2 * pi * (3 * sqrt(pi) / 4 + 41 / 45)

# The rational part of the expression is derived from the residues and winding numbers:
# 2 * Res(f, -1) + 1 * Res(f, -3) = 2 * (-2/5) + 1 * (-1/9) = -4/5 - 1/9
rational_part = Fraction(-4, 5) - Fraction(1, 9)
rational_numerator = abs(rational_part.numerator)
rational_denominator = rational_part.denominator

print("The imaginary part of the sum of the integrals is given by the equation:")
print(f"Im(I) = -2 * pi * (3 * sqrt(pi) / 4 + {rational_numerator} / {rational_denominator})")

print("\nThe numbers that appear in this final equation are:")
print(2)
print("pi")
print(3)
print("sqrt(pi)")
print(4)
print(rational_numerator)
print(rational_denominator)