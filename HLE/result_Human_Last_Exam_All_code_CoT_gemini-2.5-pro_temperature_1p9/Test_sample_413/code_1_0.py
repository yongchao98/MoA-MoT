from fractions import Fraction

# From the step-by-step derivation, we found the coefficients a, b, and c.
a = Fraction(-3, 2)
b = Fraction(-19, 16)
c = Fraction(27, 32)

# Define the polynomial f(x)
def f(x):
  return x**3 + a * x**2 + b * x + c

# The value we need to compute is f(3)
x_val = 3
result = f(x_val)

# The result is a fraction, so we will output its numerator and denominator.
numerator = result.numerator
denominator = result.denominator

# To meet the output requirement of printing the final equation,
# we format the string for f(3).
print(f"We have the function f(x) = x^3 + ({a.numerator}/{a.denominator})x^2 + ({b.numerator}/{b.denominator})x + ({c.numerator}/{c.denominator}).")
print(f"We want to compute f(3).")
# The calculation: f(3) = 3^3 - (3/2)*3^2 - (19/16)*3 + 27/32
# Step 1: Cubes and squares
val1 = 3**3
val2 = a * 3**2
val3 = b * 3
val4 = c
print(f"f(3) = {val1} + ({val2.numerator}/{val2.denominator}) + ({val3.numerator}/{val3.denominator}) + ({val4.numerator}/{val4.denominator})")

# Step 2: Combine terms
# 27 - 27/2 - 57/16 + 27/32
# = 864/32 - 432/32 - 114/32 + 27/32
term1_32 = val1 * 32
term2_32 = val2.numerator * (32 // val2.denominator)
term3_32 = val3.numerator * (32 // val3.denominator)
term4_32 = val4.numerator
print(f"f(3) = ({term1_32}/32) + ({term2_32}/32) + ({term3_32}/32) + ({term4_32}/32)")
final_numerator = term1_32 + term2_32 + term3_32 + term4_32
final_denominator = 32
print(f"f(3) = {final_numerator}/{final_denominator}")
print(f"The exact value of f(3) is {numerator}/{denominator}.")
