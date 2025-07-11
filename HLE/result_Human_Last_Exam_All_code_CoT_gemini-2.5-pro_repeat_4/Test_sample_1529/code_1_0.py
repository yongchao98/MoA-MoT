# The original curve y^2 = x^6 + 2x^3 + 4x^2 + 4x + 1 is a genus 2 curve.
# It is birationally equivalent to an elliptic curve because it has a rational point.
# A minimal integral model for the equivalent elliptic curve is y^2 = x^3 + x^2 - 2x - 1.
# We will calculate the discriminant of this minimal model.

# The general form of the cubic is x^3 + a*x^2 + b*x + c.
# For our equation, the coefficients of the polynomial on the right-hand side are:
a = 1
b = -2
c = -1

# The discriminant for a curve y^2 = x^3 + ax^2 + bx + c is given by the formula:
# Delta = a^2*b^2 - 4*b^3 - 4*a^3*c - 27*c^2 + 18*a*b*c

# Calculate each term of the discriminant formula
term1 = a**2 * b**2
term2 = -4 * b**3
term3 = -4 * a**3 * c
term4 = -27 * c**2
term5 = 18 * a * b * c

# Calculate the final discriminant value
discriminant = term1 + term2 + term3 + term4 + term5

# Print the final equation and the step-by-step calculation of its discriminant.
print(f"The minimal model of the elliptic curve is y^2 = x^3 + ({a})x^2 + ({b})x + ({c}).")
print("The discriminant is calculated as follows:")
print(f"  a^2*b^2      = ({a})^2 * ({b})^2         = {term1}")
print(f"  -4*b^3       = -4 * ({b})^3          = {term2}")
print(f"  -4*a^3*c     = -4 * ({a})^3 * ({c})      = {term3}")
print(f"  -27*c^2      = -27 * ({c})^2         = {term4}")
print(f"  +18*a*b*c    = 18 * ({a}) * ({b}) * ({c}) = {term5}")
print("--------------------------------------------------")
print(f"Minimal Discriminant = {term1} + {term2} + {term3} + {term4} + {term5} = {discriminant}")