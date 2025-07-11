# Define variables as strings for printing the expression.
# The variable 'z' represents the complex variable in the problem.
z = 'z'
# The variable 'w' represents the principal cube root of unity, exp(2*pi*i/3).
w = 'w'

# The numerator of the final expression is 8.
numerator = 8

# The denominator is composed of a polynomial part and a Gamma function part.
# The polynomial part comes from the n=1 and n=2 terms of the product.
# (1 - z**3/1**3) * (1 - z**3/2**3) = (1 - z**3) * (1 - z**3/8) = (1-z**3)*(8-z**3)/8
# The denominator of this is the numerator of the final expression.
denom_poly_part = f"(1 - {z}**3) * (8 - {z}**3)"

# The Gamma function part comes from the infinite product identity.
denom_gamma_part = f"Gamma(1 - {z}) * Gamma(1 - {w}*{z}) * Gamma(1 - {w}**2*{z})"

# Construct the complete final expression as a string.
final_expression = f"{numerator} / ({denom_poly_part} * {denom_gamma_part})"

# Print the final result with explanations for the notation.
print("The value of the infinite product is given by the expression:")
print(final_expression)
print("\nWhere:")
print(" - z is a complex variable.")
print(" - Gamma(x) denotes the Gamma function.")
print(f" - w represents the principal cube root of unity, w = exp(2*pi*i/3).")
print("\nThe numbers in the final equation are 1, 3, 8 as shown.")
