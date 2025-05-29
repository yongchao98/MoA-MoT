from fractions import Fraction

# Given fraction
numerator = 869
denominator = 1020

# Use the Fraction class to automatically reduce the fraction
reduced_fraction = Fraction(numerator, denominator)

print(f"{reduced_fraction.numerator}/{reduced_fraction.denominator}")