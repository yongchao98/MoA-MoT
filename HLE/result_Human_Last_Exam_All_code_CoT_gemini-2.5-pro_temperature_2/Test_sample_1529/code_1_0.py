import sympy

# Define the variable for the polynomial
x = sympy.Symbol('x')

# The curve is defined by y^2 = f(x). Let's define the polynomial f(x).
# The polynomial is x^6 + 2x^3 + 4x^2 + 4x + 1.
# We can represent it by its coefficients for each power of x from 6 down to 0.
c6 = 1
c5 = 0
c4 = 0
c3 = 2
c2 = 4
c1 = 4
c0 = 1

f = c6*x**6 + c5*x**5 + c4*x**4 + c3*x**3 + c2*x**2 + c1*x**1 + c0*x**0

# The "minimal discriminant" for this curve corresponds to the discriminant of the
# polynomial f(x), as the given model can be shown to be minimal.
# A minimal model is one that cannot be simplified further by a change of variables
# while keeping integer coefficients.

# We calculate the discriminant of the polynomial f(x).
discriminant_value = sympy.discriminant(f, x)

# Print the final equation with all its numbers, as requested.
print(f"The curve is defined by the equation y^2 = {c6}*x^6 + {c3}*x^3 + {c2}*x^2 + {c1}*x + {c0}")
print(f"The minimal discriminant of this curve is:")
print(discriminant_value)