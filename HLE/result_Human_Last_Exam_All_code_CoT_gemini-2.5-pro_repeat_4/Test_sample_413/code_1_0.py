from fractions import Fraction

# Step 1: Determine the coefficients a, b, and c.
# Based on the analysis, we found a = -3/4, b = -13/16.
a = Fraction(-3, 4)
b = Fraction(-13, 16)

# The condition f(k-1)f(k+1) >= 0 for all integers k implies
# that f(x) has an integer root. We found this root to be x=2.
# f(2) = 0 => 2**3 + a*(2**2) + b*2 + c = 0
# 8 + 4a + 2b + c = 0
# c = -8 - 4a - 2b
c = -8 - 4 * a - 2 * b

# The function is f(x) = x^3 + ax^2 + bx + c
print(f"The polynomial is f(x) = x^3 + ({a.numerator}/{a.denominator})x^2 + ({b.numerator}/{b.denominator})x + ({c.numerator}/{c.denominator})")

# Step 2: Compute f(3)
x = 3
f_3 = x**3 + a * x**2 + b * x + c

# The final equation to calculate is f(3) = 3^3 + a*3^2 + b*3 + c
print(f"We want to compute the value of f(3).")
print(f"f(3) = 3^3 + ({a.numerator}/{a.denominator})*(3^2) + ({b.numerator}/{b.denominator})*3 + ({c.numerator}/{c.denominator})")

# Let's show the calculation step by step with fractions
term1 = Fraction(3**3)
term2 = a * Fraction(3**2)
term3 = b * Fraction(3)
term4 = c

print(f"f(3) = {term1.numerator}/{term1.denominator} + ({term2.numerator}/{term2.denominator}) + ({term3.numerator}/{term3.denominator}) + ({term4.numerator}/{term4.denominator})")

# Combine terms
result = term1 + term2 + term3 + term4
print(f"f(3) = {result.numerator}/{result.denominator}")