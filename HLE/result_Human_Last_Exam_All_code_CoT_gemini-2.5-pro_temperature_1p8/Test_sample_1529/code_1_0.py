# Coefficients of the minimal Weierstrass model y^2 = x^3 + Ax + B
A = -7
B = 20

# Formula for the discriminant of an elliptic curve
discriminant = -16 * (4 * A**3 + 27 * B**2)

# Calculation of each term in the formula
term1 = 4 * A**3
term2 = 27 * B**2

print(f"The elliptic curve model is y^2 = x^3 + {A}x + {B}")
print(f"The discriminant is calculated as -16 * (4 * ({A})^3 + 27 * ({B})^2)")
print(f"The first term in the parenthesis is 4 * {A**3} = {term1}")
print(f"The second term in the parenthesis is 27 * {B**2} = {term2}")
print(f"The full expression is -16 * ({term1} + {term2}) = {discriminant}")