import fractions

# The user provides an integer k >= 2.
# We will calculate the value of the limit lim_{m -> infinity} ln(f(m)) / ln(m).
# Based on the reasoning in the plan, the value of this limit is 1 - 1/(2*k).
# This Python code calculates and displays this value for a given k.

# You can change this value to any integer >= 2.
k = 2

# Using the fractions module to get an exact rational number for the result.
one_int = 1
term_denominator = 2 * k
result_fraction = fractions.Fraction(term_denominator - 1, term_denominator)

print(f"The problem is to compute the limit L = lim_{m->inf} ln(f(m))/ln(m) for k = {k}.")
print("Based on extremal graph theory, the value of the limit is given by the formula 1 - 1/(2k).")
print("\nCalculation:")

# We output each number in the final equation as requested.
print(f"{one_int} - 1 / (2 * {k}) = {one_int} - 1 / {term_denominator} = {result_fraction.numerator}/{result_fraction.denominator}")

# Example for another k
k = 3
term_denominator = 2 * k
result_fraction = fractions.Fraction(term_denominator - 1, term_denominator)
print(f"\nFor k = {k}:")
print(f"{one_int} - 1 / (2 * {k}) = {one_int} - 1 / {term_denominator} = {result_fraction.numerator}/{result_fraction.denominator}")