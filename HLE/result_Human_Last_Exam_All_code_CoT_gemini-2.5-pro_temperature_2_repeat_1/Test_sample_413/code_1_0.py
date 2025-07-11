from fractions import Fraction

# From our derivation, the function is f(x) = x^3 - (3/8)x^2 - (5/8)x.
# So, a = -3/8, b = -5/8, c = 0.
a = Fraction(-3, 8)
b = Fraction(-5, 8)
c = Fraction(0)

# We want to compute the value of f(3).
x = 3

# Calculation of f(3)
term1 = x**3
term2 = a * (x**2)
term3 = b * x
result = term1 + term2 + term3 + c

# Printing the calculation steps
print(f"We compute f(3) for f(x) = x^3 + ({a})x^2 + ({b})x + {c}:")
print(f"f(3) = ({x})^3 + ({a})*({x})^2 + ({b})*({x}) + {c}")
print(f"f(3) = {term1} + {term2} + {term3} + {c}")
print(f"f(3) = {term1} + ({term2 + term3})")
term1_frac = Fraction(term1)
print(f"f(3) = {term1_frac} - {abs(term2+term3)}")
common_denom_term1 = Fraction(term1_frac.numerator * result.denominator, result.denominator)
print(f"f(3) = {common_denom_term1} - {abs(term2+term3)}")
print(f"f(3) = {result.numerator}/{result.denominator}")