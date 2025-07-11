import sympy

# Define the variable for the polynomial
x = sympy.Symbol('x')

# The curve is defined by the equation y^2 = f(x).
# The polynomial f(x) is x^6 + 2*x^3 + 4*x^2 + 4*x + 1.
# The numbers defining this polynomial are the coefficients:
# c6=1, c3=2, c2=4, c1=4, c0=1.
f = x**6 + 2*x**3 + 4*x**2 + 4*x + 1

# The minimal discriminant of the curve is the discriminant of the polynomial f(x).
# We compute this value using sympy.
minimal_discriminant = sympy.discriminant(f, x)

# The final result is an integer. The equation for the discriminant is:
# Discriminant = -1728
# We print the number from this final equation.
print(minimal_discriminant)