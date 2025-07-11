from fractions import Fraction

# From the analysis, the function is f(x) = (x-n)^3 - (x-n) with n = 1/4.
# We want to compute f(3).
n = Fraction(1, 4)
x = 3

# Define the function
def f(x_val, n_val):
    return (x_val - n_val)**3 - (x_val - n_val)

# Calculate f(3)
result = f(Fraction(x), n)

# The result is a fraction. We need to output the final equation.
# f(3) = (3 - 1/4)^3 - (3 - 1/4)
#      = (11/4)^3 - 11/4
#      = 1331/64 - 11/4
#      = 1331/64 - 176/64
#      = 1155/64
term1_val = (Fraction(x) - n)**3
term2_val = (Fraction(x) - n)
term1_base = Fraction(x) - n

print(f"The function is f(x) = (x - 1/4)^3 - (x - 1/4).")
print(f"We want to compute f(3).")
print(f"f(3) = (3 - 1/4)^3 - (3 - 1/4)")
print(f"     = ({term1_base.numerator}/{term1_base.denominator})^3 - ({term2_val.numerator}/{term2_val.denominator})")
print(f"     = {term1_val.numerator}/{term1_val.denominator} - {term2_val.numerator*16}/{term2_val.denominator*16}")
print(f"     = {(term1_val.numerator - term2_val.numerator*16)}/{term1_val.denominator}")
print(f"The exact value of f(3) is: {result.numerator}/{result.denominator}")