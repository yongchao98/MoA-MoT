from fractions import Fraction

# From the analysis, we found that alpha = 1/4.
alpha = Fraction(1, 4)

# The roots of the polynomial f(x) are alpha-1, alpha, and alpha+1.
r1 = alpha - 1
r2 = alpha
r3 = alpha + 1

# The value we want to compute is f(3).
x = 3

# Calculate the terms (x-r_i) for f(x) = (x-r1)(x-r2)(x-r3).
term1 = x - r1
term2 = x - r2
term3 = x - r3

# Compute the final result.
result = term1 * term2 * term3

# Print the final equation with all its numbers as requested.
print(f"f(3) = (3 - ({r1})) * (3 - ({r2})) * (3 - ({r3}))")
print(f"f(3) = ({term1}) * ({term2}) * ({term3})")
print(f"{term1.numerator}/{term1.denominator} * {term2.numerator}/{term2.denominator} * {term3.numerator}/{term3.denominator} = {result.numerator}/{result.denominator}")
print(f"The exact value of f(3) is {result}.")