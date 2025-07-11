import fractions

# From the derivation, the coefficients of the polynomial f(x) = x^3 + ax^2 + bx + c are:
# a = -3/8
# b = -5/8
# c = 0
# The function is f(x) = x^3 - (3/8)x^2 - (5/8)x.
# We need to compute the exact value of f(3).

# Define the coefficients as fractions for exact arithmetic
a = fractions.Fraction(-3, 8)
b = fractions.Fraction(-5, 8)
c = fractions.Fraction(0, 1)

# The value of x to evaluate
x = 3

# Calculate f(3)
# f(3) = 3^3 + a*3^2 + b*3 + c
# f(3) = 27 + a*9 + b*3 + c
f_3 = x**3 + a * x**2 + b * x + c

# The problem asks to output the final equation with each number.
# f(3) = 27 - 27/8 - 15/8 = 27 - 42/8 = 27 - 21/4 = (108 - 21)/4 = 87/4
term1 = x**3
term2 = a * x**2
term3 = b * x
term4 = c

print(f"The function is f(x) = x^3 + ({a.numerator}/{a.denominator})x^2 + ({b.numerator}/{b.denominator})x + {c.numerator}")
print(f"We want to compute f(3).")
print(f"f(3) = 3^3 + ({a.numerator}/{a.denominator})*3^2 + ({b.numerator}/{b.denominator})*3 + {c.numerator}")
print(f"f(3) = {term1} + ({term2.numerator}/{term2.denominator}) + ({term3.numerator}/{term3.denominator}) + {term4.numerator}")
print(f"f(3) = {term1} - {abs(term2.numerator)}/{term2.denominator} - {abs(term3.numerator)}/{term3.denominator}")
# Combine the fractional terms
sum_of_fractions = term2 + term3
print(f"f(3) = {term1} - {abs(sum_of_fractions.numerator)}/{sum_of_fractions.denominator}")
# Combine all terms
final_value = term1 + sum_of_fractions
print(f"f(3) = {final_value.numerator}/{final_value.denominator}")

# Final answer in the required format
print(f"<<<{f_3.numerator}/{f_3.denominator}>>>")