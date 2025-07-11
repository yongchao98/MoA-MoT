from fractions import Fraction

# From our derivation, the coefficients are a = -3/8, b = -5/8, c = 0.
a = Fraction(-3, 8)
b = Fraction(-5, 8)
c = Fraction(0, 1)

# The function is f(x) = x^3 + ax^2 + bx + c
# We want to compute f(3)
x = 3

# Calculate each term of f(3)
term1 = x**3
term2 = a * (x**2)
term3 = b * x
term4 = c

# Calculate the final value
f_3 = term1 + term2 + term3 + term4

# Print the calculation steps with exact fractions
print(f"The function is f(x) = x^3 + ({a})x^2 + ({b})x + {c}")
print(f"We want to compute f(3):")
print(f"f(3) = 3^3 + ({a}) * 3^2 + ({b}) * 3 + {c}")
print(f"f(3) = {term1} + ({a}) * 9 + ({b}) * 3 + {c}")
print(f"f(3) = {term1} + ({term2}) + ({term3}) + {term4}")
print(f"f(3) = {term1} - {abs(term2)} - {abs(term3)}")

# To combine the fractions
common_denominator = f_3.denominator
term1_frac = Fraction(term1 * common_denominator, common_denominator)
term2_frac = term2
term3_frac = term3

print(f"f(3) = {term1_frac} - {abs(term2_frac)} - {abs(term3_frac)}")
sum_of_neg_terms = abs(term2_frac) + abs(term3_frac)
print(f"f(3) = {term1_frac} - {sum_of_neg_terms}")

numerator = term1_frac.numerator - sum_of_neg_terms.numerator
print(f"f(3) = Fraction({numerator}, {common_denominator})")

# Print the final result
print(f"The exact value of f(3) is: {f_3}")