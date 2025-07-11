from fractions import Fraction

# From our derivation, we found the function must be of the form f(x) = x(x-1)(x-beta)
# which corresponds to a = -(1+beta) and b = beta.
# Using the derivative condition b - a/2 = -7/16, we solved for beta.
beta = Fraction(-5, 8)
a = -(1 + beta)
b = beta
c = 0

print(f"The coefficients of the polynomial f(x) = x^3 + ax^2 + bx + c are:")
print(f"a = {a}")
print(f"b = {b}")
print(f"c = {c}")
print(f"So, f(x) = x^3 - (3/8)x^2 - (5/8)x")

# We are asked to compute the value of f(3).
x_val = 3
f_val = Fraction(x_val)**3 + a * Fraction(x_val)**2 + b * Fraction(x_val)

print("\nNow we compute f(3):")
term1 = Fraction(x_val**3)
term2 = a * Fraction(x_val**2)
term3 = b * Fraction(x_val)

print(f"f(3) = 3^3 - (3/8)*3^2 - (5/8)*3")
print(f"f(3) = {term1} - {abs(term2)} - {abs(term3)}")

# Show the combination of fractional terms
term23_sum = term2 + term3
print(f"f(3) = {term1} - ({abs(term2)} + {abs(term3)}) = {term1} - {abs(term23_sum)}")

# Final calculation in a common denominator
term1_common_denom = term1.limit_denominator(f_val.denominator)

print(f"The final equation is: {term1_common_denom.numerator}/{term1_common_denom.denominator} - {abs(term23_sum.numerator)}/{term23_sum.denominator} = {f_val.numerator}/{f_val.denominator}")
print(f"So the final answer is: {f_val.numerator}/{f_val.denominator}")