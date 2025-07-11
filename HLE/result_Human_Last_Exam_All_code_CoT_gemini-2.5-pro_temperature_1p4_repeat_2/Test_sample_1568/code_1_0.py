import math

# The problem is to evaluate the infinite product:
# Product_{n=3 to infinity} (1 - z^3 / n^3)
#
# The final result can be expressed using the Gamma function.
# The following code prints the symbolic formula for the product.

# The equation is of the form: Product = Numerator / Denominator
# Let's define the parts of the formula.
# In the formula, 'z' is a complex variable.
# 'Gamma' represents the Gamma function.
# 'exp(j*theta)' is the complex exponential e^(j*theta).

# The formula contains several numbers: 1, 3, 8, 2, 4.
# As requested, the code will output these numbers directly in the final equation.

numerator = "1"

term1 = "(1 - z**3)"
term2 = "(1 - z**3 / 8)"
term3 = "Gamma(1 - z)"
term4 = "Gamma(1 - z * exp(i*2*pi/3))"
term5 = "Gamma(1 - z * exp(i*4*pi/3))"

denominator = f"{term1} * {term2} * {term3} * {term4} * {term5}"

print("The infinite product is equal to:")
print(f"{numerator} / ({denominator})")
print("\nWhere:")
print("'z' is a complex variable.")
print("'Gamma' is the Gamma function.")
print("'i' is the imaginary unit.")
print("'pi' is the mathematical constant pi.")
print("'exp' is the exponential function.")
