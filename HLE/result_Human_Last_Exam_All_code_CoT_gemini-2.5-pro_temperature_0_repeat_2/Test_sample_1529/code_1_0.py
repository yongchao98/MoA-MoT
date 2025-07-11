import sympy

# The problem is to find the minimal discriminant of the curve y^2 = x^6 + 2*x^3 + 4*x^2 + 4*x + 1.
# We found that the curve is birationally equivalent to y^2 = g(x), where g(x) is a simpler polynomial.
# The final equation for the minimal model is y^2 = x^5 - 6*x^4 + 15*x^3 - 18*x^2 + 13*x - 4.
# We will now compute the discriminant of this polynomial.

# Define the symbolic variable
x = sympy.Symbol('x')

# The polynomial from the minimal model
# The numbers in this equation are: 1, -6, 15, -18, 13, -4
g = x**5 - 6*x**4 + 15*x**3 - 18*x**2 + 13*x - 4

# Calculate the discriminant of the polynomial g(x)
minimal_discriminant = sympy.discriminant(g, x)

# Print the result
print(minimal_discriminant)