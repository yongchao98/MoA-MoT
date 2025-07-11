from fractions import Fraction

# Coefficients determined from the problem conditions
a = Fraction(-3, 4)
b = Fraction(-13, 16)
c = Fraction(0, 1) # This is c=0 based on the analysis

# The function f(x)
def f(x):
  x = Fraction(x)
  return x**3 + a * x**2 + b * x + c

# Value to compute
x_val = 3

# Calculation steps
val_f3 = f(x_val)
term1 = x_val**3
term2 = a * x_val**2
term3 = b * x_val

print(f"Based on the problem analysis, the function is f(x) = x^3 - (3/4)x^2 - (13/16)x.")
print(f"We need to compute f(3).")
print(f"f(3) = 3^3 + ({a}) * 3^2 + ({b}) * 3 + {c}")
print(f"f(3) = {term1} + ({a}) * 9 + ({b}) * 3")
print(f"f(3) = {term1} + ({term2}) + ({term3})")
print(f"f(3) = {term1} - {abs(term2)} - {abs(term3)}")

# Convert to common denominator for display
term1_common = term1.limit_denominator(16)
term2_common = term2.limit_denominator(16)
term3_common = term3.limit_denominator(16)

print(f"f(3) = {term1_common.numerator}/{term1_common.denominator} - {abs(term2_common.numerator)}/{abs(term2_common.denominator)} - {abs(term3_common.numerator)}/{abs(term3_common.denominator)}")

sum_numerator = term1_common.numerator + term2_common.numerator + term3_common.numerator

print(f"f(3) = ({term1_common.numerator} - {abs(term2_common.numerator)} - {abs(term3_common.numerator)}) / 16")
print(f"f(3) = {sum_numerator}/16")
print(f"The exact value is {val_f3}")
<<<285/16>>>